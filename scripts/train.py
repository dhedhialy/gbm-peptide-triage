#!/usr/bin/env python3
"""Train the two sequence models (random forest on transparent descriptors),
calibrate them with split-conformal prediction, benchmark random vs
sequence-similarity splits, and export everything the web app needs.

Outputs:
  public/model.json     - forest structure, importances, calibration, fixtures
  public/benchmark.json - protocol metrics, reliability, alpha sweeps

Run: python3 scripts/train.py
"""
import json
import math
import random
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "public"

ALPHA = 0.10          # conformal: target 1-alpha marginal coverage
SIM_CLUSTER = 0.80    # sequences with sim >= this join one cluster
OOD_SIM_FLOOR = 0.40  # floor for the OOD gate
OOD_SIM_PCT = 0.10    # the gate is the p10 of held-out similarity, floored above
SPLIT = (0.60, 0.20, 0.20)  # train / calibrate / test, by cluster
SEEDS = [0, 1, 2, 3, 4]
SWEEP_ALPHAS = [0.40, 0.30, 0.20, 0.15, 0.10, 0.05]  # exported so the UI can move the abstention dial
RF_PARAMS = dict(n_estimators=200, max_depth=6, min_samples_leaf=3,
                 class_weight="balanced", random_state=0, n_jobs=-1)
N_TREES_DEPLOY = 200

# ---------------------------------------------------------------- features --

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
    "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
    "Y": -1.3, "V": 4.2,
}
AA_MW = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "Q": 128.1307, "E": 129.1155, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
PKA_ACID = {"D": 3.65, "E": 4.25, "C": 8.33, "Y": 10.1}
PKA_BASE = {"K": 10.5, "R": 12.0, "H": 6.0}
PKA_NTERM, PKA_CTERM, PH = 9.6, 2.35, 7.4

BASE_FEATURES = [
    "length", "net_charge", "gravy", "aromatic_frac", "polar_frac",
    "charged_frac", "aliphatic_frac", "gly_frac", "pro_frac", "cys_frac",
    "acidic_frac", "basic_frac", "max_hydro_run", "mw",
]
FEATURES = BASE_FEATURES + [f"frac_{a}" for a in AA]
HYDRO_POS = set(a for a, v in AA_KD.items() if v > 0)


def features(seq: str) -> list[float]:
    n = len(seq)
    charge = 1.0 / (1.0 + 10 ** (PH - PKA_NTERM))
    charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - PH))
    for aa, pk in PKA_ACID.items():
        charge -= seq.count(aa) / (1.0 + 10 ** (pk - PH))
    for aa, pk in PKA_BASE.items():
        charge += seq.count(aa) / (1.0 + 10 ** (PH - pk))
    run = best = 0
    for c in seq:
        run = run + 1 if c in HYDRO_POS else 0
        best = max(best, run)
    f = [
        float(n),
        charge,
        sum(AA_KD[c] for c in seq) / n,
        sum(seq.count(a) for a in "FWY") / n,
        sum(seq.count(a) for a in "STNQ") / n,
        sum(seq.count(a) for a in "DEKRH") / n,
        sum(seq.count(a) for a in "VLI") / n,
        seq.count("G") / n,
        seq.count("P") / n,
        seq.count("C") / n,
        sum(seq.count(a) for a in "DE") / n,
        sum(seq.count(a) for a in "KRH") / n,
        float(best),
        sum(AA_MW[c] for c in seq) + 18.0153,
    ]
    f.extend(seq.count(a) / n for a in AA)
    assert len(f) == len(FEATURES)
    return f


def read_fasta(path: Path) -> list[str]:
    seqs, s = [], ""
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if s:
                seqs.append(s)
            s = ""
        else:
            s += line
    if s:
        seqs.append(s)
    return seqs


# ------------------------------------------------------------- similarity ---

_lcs_cache: dict[tuple[str, str], int] = {}


def lcs(a: str, b: str) -> int:
    key = (a, b) if a <= b else (b, a)
    if key in _lcs_cache:
        return _lcs_cache[key]
    if len(a) < len(b):
        a, b = b, a
    prev = [0] * (len(b) + 1)
    for ca in a:
        cur = [0]
        for j, cb in enumerate(b, 1):
            cur.append(prev[j - 1] + 1 if ca == cb else max(prev[j], cur[j - 1]))
        prev = cur
    v = prev[-1]
    _lcs_cache[key] = v
    return v


def sim(a: str, b: str, floor: float = SIM_CLUSTER) -> float:
    """2 * LCS / (len_a + len_b) -- identical definition in src/lib/identity.ts

    `floor` enables a cheap reject: Dice-on-LCS is at most 2r/(1+r) with
    r = min/max length, so pairs whose upper bound is below the caller's
    threshold skip the LCS entirely.
    """
    if a == b:
        return 1.0
    lo, hi = (len(a), len(b)) if len(a) <= len(b) else (len(b), len(a))
    if 2.0 * lo / (lo + hi) < floor:
        return 0.0
    return 2.0 * lcs(a, b) / (len(a) + len(b))


_cluster_cache: dict[tuple, list[int]] = {}


def cluster(seqs: list[str], threshold: float) -> list[int]:
    key = (threshold, tuple(seqs))
    if key in _cluster_cache:
        return _cluster_cache[key]
    parent = list(range(len(seqs)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            if sim(seqs[i], seqs[j], threshold) >= threshold:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[rj] = ri
    roots = {}
    out = [roots.setdefault(find(i), len(roots)) for i in range(len(seqs))]
    _cluster_cache[key] = out
    return out


def max_sim_to(query: str, pool: list[str]) -> float:
    # floor 0.0: exact similarity (the browser reports exact values too);
    # the OOD gate (derived in run_protocol) is applied to this raw maximum.
    return max((sim(query, p, 0.0) for p in pool), default=0.0)


# ----------------------------------------------------------- conformal ------

def conformal_quantile(n: int, alpha: float) -> float:
    """Finite-sample split-conformal quantile level (ceiling rule)."""
    return min(1.0, math.ceil((n + 1) * (1 - alpha)) / n)


def calibrate(p_true: np.ndarray, alpha: float) -> float:
    """Quantile of nonconformity 1 - p(true class) on the calibration set."""
    nc = 1.0 - np.asarray(p_true, dtype=float)
    return float(np.quantile(nc, conformal_quantile(len(nc), alpha), method="higher"))


def pred_sets(p: np.ndarray, q: float) -> np.ndarray:
    """LAC prediction sets as (n, 2) bool: both True = abstain."""
    return (1.0 - p) <= (q + 1e-12)


def set_stats(p: np.ndarray, y: np.ndarray, q: float) -> dict:
    sets = pred_sets(p, q)
    sizes = sets.sum(axis=1)
    covered = sets[np.arange(len(y)), y.astype(int)]
    singleton = sizes == 1
    empty = sizes == 0
    singleton_acc = (
        float((sets[singleton].argmax(axis=1) == y[singleton]).mean())
        if singleton.sum() else float("nan")
    )
    return {
        "coverage": float(covered.mean()),
        "abstain_rate": float((sizes == 2).mean()),
        "empty_rate": float(empty.mean()),
        "singleton_acc": singleton_acc,
        "mean_set_size": float(sizes.mean()),
    }


# ------------------------------------------------------------ protocols -----

def split_random(n, y, seed, fracs=SPLIT):
    rng = random.Random(seed)
    idx = list(range(n))
    rng.shuffle(idx)
    by_y = {0: [i for i in idx if y[i] == 0], 1: [i for i in idx if y[i] == 1]}
    out = []
    for part in by_y.values():
        n1 = int(round(fracs[0] * len(part)))
        n2 = int(round((fracs[0] + fracs[1]) * len(part)))
        out.append((part[:n1], part[n1:n2], part[n2:]))
    tr = out[0][0] + out[1][0]
    ca = out[0][1] + out[1][1]
    te = out[0][2] + out[1][2]
    rng.shuffle(tr)
    rng.shuffle(ca)
    rng.shuffle(te)
    return tr, ca, te


def split_similarity(seqs, y, fracs=SPLIT, seed=0):
    """Cluster at SIM_CLUSTER, then fill 60/20/20 inside each label bucket, so
    no test cluster shares >=80% similarity with any training cluster while
    both classes stay represented everywhere. Cluster membership never varies;
    only the seed-controlled shuffle of cluster order does."""
    labels = cluster(seqs, SIM_CLUSTER)
    groups: dict[int, dict[int, list[int]]] = {}
    for i, l in enumerate(labels):
        groups.setdefault(l, {0: [], 1: []})[y[i]].append(i)
    buckets = {0: [], 1: []}
    for members in groups.values():
        maj = 1 if len(members[1]) > len(members[0]) else 0
        buckets[maj].append(sorted(members[0] + members[1]))
    rng = random.Random(seed)
    tr, ca, te = [], [], []
    for maj in (0, 1):
        cl = buckets[maj]
        rng.shuffle(cl)
        size = sum(len(c) for c in cl)
        t_end, c_end = fracs[0] * size, (fracs[0] + fracs[1]) * size
        acc = 0
        for c in cl:
            if acc < t_end:
                tr += c
            elif acc < c_end:
                ca += c
            else:
                te += c
            acc += len(c)
    rng.shuffle(tr)
    rng.shuffle(ca)
    rng.shuffle(te)
    assert len(set(tr) | set(ca) | set(te)) == len(seqs)
    return tr, ca, te, labels


def fit(X, y, tr):
    m = RandomForestClassifier(**RF_PARAMS).fit([X[i] for i in tr], y[tr])
    return m


def run_protocol(seqs, y, protocol, seed=0):
    X = [features(s) for s in seqs]
    y = np.array(y)
    if protocol == "random":
        tr, ca, te = split_random(len(seqs), y, seed)
    else:
        tr, ca, te, labels = split_similarity(seqs, y, seed=seed)
    model = fit(X, y, tr)
    p_ca = model.predict_proba([X[i] for i in ca])
    q = calibrate(p_ca[np.arange(len(ca)), y[ca]], ALPHA)
    p_te = model.predict_proba([X[i] for i in te])
    p_te1 = p_te[:, 1]
    res = {
        "auroc": float(roc_auc_score(y[te], p_te1)),
        "auprc": float(average_precision_score(y[te], p_te1)),
        "n_train": len(tr),
        "n_cal": len(ca),
        "n_test": len(te),
        "conformal": set_stats(p_te, y[te], q),
        "q": q,
    }
    if protocol == "similarity":
        res["n_clusters"] = len(set(labels))
    # reliability diagram (10 fixed bins, on test)
    bins = np.clip((p_te1 * 10).astype(int), 0, 9)
    rel = []
    for b in range(10):
        m = bins == b
        if m.sum():
            rel.append({
                "predicted": round((b + 0.5) / 10, 3),
                "observed": round(float(y[te][m].mean()), 4),
                "count": int(m.sum()),
            })
    res["reliability"] = rel
    # alpha sweep: what abstention costs / buys
    sweep = []
    for a in SWEEP_ALPHAS:
        qa = calibrate(p_ca[np.arange(len(ca)), y[ca]], a)
        sets = pred_sets(p_te, qa)
        sizes = sets.sum(axis=1)
        conf = sizes == 1
        acc = (float((sets[conf].argmax(axis=1) == y[te][conf]).mean())
               if conf.sum() else float("nan"))
        sweep.append({
            "alpha": a,
            "abstain_rate": round(float((sizes != 1).mean()), 4),
            "error_rate": round(1.0 - acc, 4),
            "coverage": round(float(sets[np.arange(len(te)), y[te].astype(int)].mean()), 4),
        })
    res["alpha_sweep"] = sweep
    # how similar are test peptides to the training pool (OOD behaviour)
    pool = [seqs[i] for i in tr]
    idents = sorted(max_sim_to(seqs[i], pool) for i in te)
    # the OOD gate is the bottom decile of held-out similarity to the training
    # pool, floored at OOD_SIM_FLOOR. Deriving it from data beats guessing a
    # constant: a fixed 0.40 almost never fires on real peptides.
    ood_sim = round(max(OOD_SIM_FLOOR, float(idents[int(OOD_SIM_PCT * (len(idents) - 1))])), 2)
    res["ood_sim"] = ood_sim
    res["test_to_train_sim"] = {
        "p10": round(float(idents[int(0.10 * (len(idents) - 1))]), 4),
        "median": round(float(idents[len(idents) // 2]), 4),
        "p90": round(float(idents[int(0.90 * (len(idents) - 1))]), 4),
        "frac_below_ood": round(float(np.mean(np.array(idents) < ood_sim)), 4),
    }
    return res


def stats_over(seqs, y, protocol, seeds=SEEDS):
    runs = [run_protocol(seqs, y, protocol, seed=s) for s in seeds]
    aurocs = [r["auroc"] for r in runs]
    auprcs = [r["auprc"] for r in runs]
    out = {
        "auroc_mean": round(float(np.mean(aurocs)), 4),
        "auroc_sd": round(float(np.std(aurocs)), 4),
        "auprc_mean": round(float(np.mean(auprcs)), 4),
        "auprc_sd": round(float(np.std(auprcs)), 4),
        "n_seeds": len(seeds),
        "conformal": runs[0]["conformal"],   # seed-0 split, reported in detail
        "q": round(runs[0]["q"], 5),
        "reliability": runs[0]["reliability"],
        "alpha_sweep": runs[0]["alpha_sweep"],
        "test_to_train_sim": runs[0]["test_to_train_sim"],
        "ood_sim": runs[0]["ood_sim"],
        "n_train": runs[0]["n_train"],
        "n_cal": runs[0]["n_cal"],
        "n_test": runs[0]["n_test"],
    }
    if "n_clusters" in runs[0]:
        out["n_clusters"] = runs[0]["n_clusters"]
    return out


# ------------------------------------------- B3Pred's own split + audit -------

# External-validation AUROC of B3Pred's RF, quoted from Kumar et al. 2021,
# Pharmaceutics 13(8):1237, Tables 1-3 (random 80:20 splits).
B3PRED_PUBLISHED_RF_VAL_AUROC = {"d1": 0.88, "d2": 0.88, "d3": 0.90}


def b3pred_protocol():
    """Train our forest on B3Pred's published train partitions and test on
    their published valid partitions (same peptides, same split), then audit
    how much of their validation set is >=80% similar to their training set."""
    pos_tr = read_fasta(RAW / "pos_train.fasta")
    pos_va = read_fasta(RAW / "pos_valid.fasta")
    out = {}
    for d in ("d1", "d2", "d3"):
        neg_tr = read_fasta(RAW / f"{d.upper()}_neg_train.fasta")
        neg_va = read_fasta(RAW / f"{d.upper()}_neg_valid.fasta")
        # drop any cross-file duplicates so each sequence keeps one label
        pos_set = set(pos_tr) | set(pos_va)
        neg_tr = [s for s in neg_tr if s not in pos_set]
        neg_va = [s for s in neg_va if s not in pos_set]
        tr_seqs, tr_y = pos_tr + neg_tr, [1] * len(pos_tr) + [0] * len(neg_tr)
        va_seqs, va_y = pos_va + neg_va, [1] * len(pos_va) + [0] * len(neg_va)
        X_tr = [features(s) for s in tr_seqs]
        X_va = [features(s) for s in va_seqs]
        model = RandomForestClassifier(**RF_PARAMS).fit(X_tr, np.array(tr_y))
        p = model.predict_proba(X_va)[:, 1]
        sims = np.array([max_sim_to(s, tr_seqs) for s in va_seqs])
        novel = sims < SIM_CLUSTER
        res = {
            "published_rf_val_auroc": B3PRED_PUBLISHED_RF_VAL_AUROC[d],
            "our_auroc": round(float(roc_auc_score(va_y, p)), 4),
            "our_auprc": round(float(average_precision_score(va_y, p)), 4),
            "n_train": len(tr_seqs),
            "n_valid": len(va_seqs),
            "frac_valid_ge_0.8_similar_to_train": round(float((sims >= SIM_CLUSTER).mean()), 4),
            "frac_valid_ge_0.6_similar_to_train": round(float((sims >= 0.6).mean()), 4),
            "median_max_sim": round(float(np.median(sims)), 4),
        }
        if novel.sum() and len(set(np.array(va_y)[novel])) == 2:
            res["our_auroc_novel_only"] = round(
                float(roc_auc_score(np.array(va_y)[novel], p[novel])), 4)
            res["n_novel_only"] = int(novel.sum())
        else:
            res["our_auroc_novel_only"] = None
            res["n_novel_only"] = int(novel.sum())
        out[d] = res
        print(f"  [{d}] our AUROC {res['our_auroc']}  "
              f"published RF {res['published_rf_val_auroc']}  "
              f"valid>=0.8sim {res['frac_valid_ge_0.8_similar_to_train']:.1%}  "
              f"novel-only {res['our_auroc_novel_only']}")
    return out


# -------------------------------------------------------------- export ------

def export_forest(model) -> dict:
    """Compact flat-array encoding of the forest, traversed by src/lib/model.ts.

    Per tree, nodes are stored in sklearn's depth-first order:
      left[i], right[i]  - node indices (-1 for leaf)
      feature[i], thr[i] - split feature index and threshold
      prob[i]            - P(class=1) at this node (leaf prediction)
    """
    trees = []
    for est in model.estimators_:
        t = est.tree_
        counts = t.value[:, 0, :]  # (n_nodes, n_classes) weighted class counts
        trees.append({
            "left": [int(v) for v in t.children_left],
            "right": [int(v) for v in t.children_right],
            "feature": [int(v) for v in t.feature],
            "thr": [float(v) for v in t.threshold],
            "prob": [float(v) for v in counts[:, 1] / counts.sum(axis=1)],
        })
    return {"n_trees": len(trees), "trees": trees,
            "importance": [round(float(v), 5) for v in model.feature_importances_]}


def prob_of(forest: dict, x: list[float]) -> float:
    """Reference forest inference (mirrors the TS implementation)."""
    total = 0.0
    for t in forest["trees"]:
        i = 0
        while t["left"][i] != -1:
            i = t["left"][i] if x[t["feature"][i]] <= t["thr"][i] else t["right"][i]
        total += t["prob"][i]
    return total / forest["n_trees"]


def local_attribution(forest: dict, x: list[float]) -> list[list]:
    """Per-path probability deltas: [[feature_index, delta], ...] sorted by |delta|.

    Walks each tree from root to leaf, recording how much each split moves the
    node probability; averages across trees (same mechanics as src/lib/model.ts).
    """
    acc: dict[int, float] = {}
    for t in forest["trees"]:
        i = 0
        while t["left"][i] != -1:
            j = t["left"][i] if x[t["feature"][i]] <= t["thr"][i] else t["right"][i]
            feat = t["feature"][i]
            acc[feat] = acc.get(feat, 0.0) + (t["prob"][j] - t["prob"][i])
            i = j
    items = [[f, d / forest["n_trees"]] for f, d in acc.items() if abs(d) > 1e-9]
    items.sort(key=lambda kv: -abs(kv[1]))
    return items


# ------------------------------------------------------------------ main ----

def main():
    random.seed(0)
    np.random.seed(0)

    pos = sorted(set(read_fasta(RAW / "Positive_B3PPs.fasta")))
    cpp = sorted(set(read_fasta(RAW / "Negative_CPPs.fasta")))
    rnd = sorted(set(read_fasta(RAW / "Negative_Random_Balanced.fasta")))

    # one sequence, one label per task (files overlap slightly)
    pos_set, cpp_set = set(pos), set(cpp)
    cpp = [s for s in cpp if s not in pos_set]
    rnd = [s for s in rnd if s not in pos_set and s not in cpp_set]

    tasks_data = {
        "bbb": (pos + cpp + rnd, [1] * len(pos) + [0] * (len(cpp) + len(rnd))),
        "cpp": (cpp + rnd, [1] * len(cpp) + [0] * len(rnd)),
    }

    benchmark = {
        "generated": "scripts/train.py",
        "model": {k.replace("__", ""): v for k, v in RF_PARAMS.items() if k != "n_jobs"},
        "alpha": ALPHA,
        "similarity": {"definition": "2*LCS/(len_a+len_b)",
                       "cluster_threshold": SIM_CLUSTER, "ood_threshold": None},
        "dataset": {"bbb_pos": len(pos), "cpp_pos": len(cpp), "random_neg": len(rnd)},
        "tasks": {},
    }
    model_json = {
        "version": 2,
        "features": FEATURES,
        "alpha": ALPHA,
        "ood_sim": None,
        "tasks": {},
        "fixtures": [],
    }

    ood_by_task = {}
    for name, (seqs, y) in tasks_data.items():
        print(f"[{name}] n={len(seqs)} pos={sum(y)}")
        for protocol in ("random", "similarity"):
            st = stats_over(seqs, y, protocol)
            print(f"  {protocol:10s} AUROC {st['auroc_mean']:.3f}±{st['auroc_sd']:.3f}"
                  f"  AUPRC {st['auprc_mean']:.3f}")
            if protocol == "similarity":
                print(f"  conformal: {st['conformal']}")
                print(f"  OOD gate: {st['ood_sim']} (p10 of held-out similarity)")
                ood_by_task[name] = st["ood_sim"]
            benchmark["tasks"].setdefault(name, {})[protocol] = st

        # deploy: similarity split, seed 0 (test clusters never trained on)
        X = [features(s) for s in seqs]
        yy = np.array(y)
        tr, ca, te, labels = split_similarity(seqs, yy, seed=0)
        model = fit(X, yy, tr)
        forest = export_forest(model)
        p_ca = model.predict_proba([X[i] for i in ca])
        nc_ca = p_ca[np.arange(len(ca)), yy[ca]]
        q = calibrate(nc_ca, ALPHA)
        qs = {str(a): calibrate(nc_ca, a) for a in SWEEP_ALPHAS}
        p_all = model.predict_proba(X)

        model_json["tasks"][name] = {
            **forest,
            "q": float(q),
            "qs": qs,
            "train_seqs": [{"s": seqs[i], "y": int(yy[i])} for i in tr],
            "train_len_range": [int(min(len(seqs[i]) for i in tr)),
                                int(max(len(seqs[i]) for i in tr))],
            "n_cal": len(ca),
            "n_test": len(te),
            "test_indices": [int(i) for i in te],
        }
        # fixtures: first 3 positives + 2 negatives of the training part
        want = [i for i in tr if yy[i] == 1][:3] + [i for i in tr if yy[i] == 0][:2]
        for i in want:
            # self-check: exported forest must reproduce sklearn exactly
            assert abs(prob_of(forest, X[i]) - p_all[i][1]) < 1e-9, "forest export mismatch"
            sets = pred_sets(p_all[i].reshape(1, -1), q)[0]
            model_json["fixtures"].append({
                "task": name,
                "s": seqs[i],
                "features": [round(v, 8) for v in X[i]],
                "prob": round(float(p_all[i][1]), 8),
                "set": [bool(sets[0]), bool(sets[1])],
                "sim_to_train": round(max_sim_to(seqs[i], [seqs[j] for j in tr]), 4),
                "top_attr": [[f, round(d, 6)] for f, d in local_attribution(forest, X[i])[:3]],
            })
        print(f"  deploy: train={len(tr)} cal={len(ca)} test={len(te)} q={q:.4f}"
              f" trees={forest['n_trees']}")

    # one OOD gate for the whole UI, derived from the primary (BBB) task
    ood_sim = ood_by_task["bbb"]
    model_json["ood_sim"] = ood_sim
    benchmark["similarity"]["ood_threshold"] = ood_sim
    print(f"OOD gate deployed: {ood_sim}")

    print("[b3pred published splits]")
    benchmark["b3pred_published_splits"] = {
        "citation": "Kumar V. et al., Pharmaceutics 2021, 13(8):1237, Tables 1-3",
        "note": ("RF external-validation AUROC as published (random 80:20 splits). "
                 "We retrain our model on their exact train/valid partitions, then "
                 "audit similarity between their valid and train peptides."),
        "datasets": b3pred_protocol(),
    }

    OUT.mkdir(exist_ok=True)
    (OUT / "model.json").write_text(json.dumps(model_json))
    (OUT / "benchmark.json").write_text(json.dumps(benchmark, indent=1))
    print(f"wrote public/model.json ({(OUT / 'model.json').stat().st_size // 1024} KB)"
          f" and public/benchmark.json")


if __name__ == "__main__":
    main()
