import { groupedBars, lineChart, type BarGroup } from "../viz/charts";
import type { AppState } from "./context";
import { esc, pct } from "./dom";

export function renderValidation(host: HTMLElement, state: AppState): void {
  const b = state.benchmark;
  if (!b) {
    host.innerHTML = `<div class="empty"><p>Loading benchmark…</p></div>`;
    return;
  }
  const d = b.b3pred_published_splits.datasets;
  const names: Record<string, string> = {
    d1: "D1 · B3PP vs CPP",
    d2: "D2 · vs balanced random",
    d3: "D3 · vs full random",
  };
  const groupsA: BarGroup[] = Object.keys(d).map((k) => ({
    label: names[k] ?? k,
    values: [
      d[k].published_rf_val_auroc,
      d[k].our_auroc,
      d[k].our_auroc_novel_only ?? 0,
    ],
  }));
  const annotations = Object.keys(d).map(
    (k) => `${pct(d[k]["frac_valid_ge_0.8_similar_to_train"])} ≥80% like train`,
  );

  const groupsB: BarGroup[] = [
    { label: "BBB · random", values: [b.tasks.bbb.random.auroc_mean, b.tasks.bbb.random.auprc_mean] },
    { label: "BBB · similar", values: [b.tasks.bbb.similarity.auroc_mean, b.tasks.bbb.similarity.auprc_mean] },
    { label: "Cell · random", values: [b.tasks.cpp.random.auroc_mean, b.tasks.cpp.random.auprc_mean] },
    { label: "Cell · similar", values: [b.tasks.cpp.similarity.auroc_mean, b.tasks.cpp.similarity.auprc_mean] },
  ];

  const sweep = (key: "bbb" | "cpp", proto: "random" | "similarity") =>
    b.tasks[key][proto].alpha_sweep
      .filter((s) => Number.isFinite(s.error_rate))
      .map((s) => ({ x: s.abstain_rate, y: s.error_rate }));

  host.innerHTML = `
    <section class="card">
      <h2>Does it generalise? — the split matters more than the model</h2>
      <p class="lede">B3Pred is the closest published predictor on this exact data. We trained
      our model on <em>their</em> published training partitions and tested on <em>their</em>
      validation partitions, so the test peptides are identical. Then we measured how many of
      those validation peptides are already near-duplicates of training peptides
      (≥80% sequence similarity).</p>
      <div id="bars-a"></div>
      <p class="hint">Note how D3 (full random decoys, 10:1 imbalance) keeps AUROC 0.92 but
      its AUPRC drops to ${d.d3.our_auprc.toFixed(2)} — AUROC hides class imbalance, which is
      exactly why the abstention accounting above is the more useful number.</p>
      <p class="hint">${esc(b.b3pred_published_splits.citation)} — published RF AUROC is
      quoted from their paper (their random 80:20 split); “ours” and “ours, novel-only” are
      measured here. Novel-only = validation peptides below 80% similarity to every training
      peptide (n = ${Object.keys(d)
        .map((k) => d[k].n_novel_only)
        .join(" / ")}).</p>
      <table class="mini-table wide">
        <thead><tr><th>split</th><th>published RF</th><th>ours (same split)</th>
        <th>ours, novel-only</th><th>valid ≥80% like train</th><th>valid ≥60% like train</th></tr></thead>
        <tbody>${Object.keys(d)
          .map(
            (k) => `<tr><td>${esc(names[k] ?? k)}</td>
          <td class="num">${d[k].published_rf_val_auroc.toFixed(2)}</td>
          <td class="num">${d[k].our_auroc.toFixed(3)}</td>
          <td class="num">${d[k].our_auroc_novel_only?.toFixed(3) ?? "—"}</td>
          <td class="num">${pct(d[k]["frac_valid_ge_0.8_similar_to_train"], 1)}</td>
          <td class="num">${pct(d[k]["frac_valid_ge_0.6_similar_to_train"], 1)}</td></tr>`,
          )
          .join("")}</tbody>
      </table>
    </section>

    <section class="card">
      <h2>Our own protocol: random split vs 80%-similarity cluster split</h2>
      <p class="lede">Random 80:20 splits leak near-duplicates between train and test, as the
      audit above shows. Clustering every peptide at 80% similarity and holding out whole
      clusters gives an honest estimate. Mean ± sd over ${b.tasks.bbb.random.n_seeds} seeds.</p>
      <div id="bars-b"></div>
      <p class="hint">BBB similarity split: ${b.tasks.bbb.similarity.n_clusters} clusters;
      test-set similarity to train pool — median ${pct(b.tasks.bbb.similarity.test_to_train_sim.median)},
      ${pct(b.tasks.bbb.similarity.test_to_train_sim.frac_below_ood, 1)} of test peptides below the
      ${pct(b.similarity.ood_threshold)} out-of-domain flag.</p>
    </section>

    <section class="card">
      <h2>The abstention dial, measured</h2>
      <p class="lede">Every point is a different coverage target α, computed on the held-out
      test clusters. Abstaining is not failure — it is the knob that buys a bounded error rate
      on the answers that remain.</p>
      <div class="grid-2">
        <div id="line-abstain"></div>
        <div id="line-reliability"></div>
      </div>
      <p class="hint">Left: error rate among <em>answered</em> peptides vs how many are left
      unanswered. The curve is not monotone step to step — with ~150 calibration peptides each
      α lands on a different empirical quantile, and the middle of the grid is noisy. The ends
      are the reliable statement: refusing ~60% of peptides instead of ~30% roughly halves the
      error rate. Right: reliability diagram (predicted probability vs observed frequency,
      diagonal = perfect calibration); the ~5-bin histogram is coarse for the same
      small-sample reason.</p>
    </section>

    <section class="card">
      <h2>Conformal accounting</h2>
      <table class="mini-table wide">
        <thead><tr><th>task</th><th>split</th><th>coverage (target ${pct(1 - b.alpha)})</th>
        <th>abstained</th><th>accuracy when answering</th><th>test n</th></tr></thead>
        <tbody>
        ${(["bbb", "cpp"] as const)
          .flatMap((t) =>
            (["random", "similarity"] as const).map((p) => {
              const c = b.tasks[t][p].conformal;
              return `<tr><td>${t === "bbb" ? "BBB penetration" : "Cell penetration"}</td>
              <td>${p}</td><td class="num">${pct(c.coverage, 1)}</td>
              <td class="num">${pct(c.abstain_rate, 1)}</td>
              <td class="num">${Number.isFinite(c.singleton_acc) ? pct(c.singleton_acc, 1) : "—"}</td>
              <td class="num">${b.tasks[t][p].n_test}</td></tr>`;
            }),
          )
          .join("")}
        </tbody>
      </table>
      <p class="hint">Coverage is empirical on the test clusters; the split-conformal
      guarantee holds marginally on exchangeable data (Lei et al., 2018), and is explicitly
      voided for peptides flagged out-of-domain.</p>
    </section>

    <section class="card warn">
      <h2>What this benchmark does not prove</h2>
      <ul class="reasons">
        <li>We could not re-run B3Pred itself under similarity splits: its standalone package
        references a server-side model file (<code>/gpsr/models/bp3red/RF_model</code>) that is
        not distributed. Published numbers are quoted, not reproduced.</li>
        <li>The BBB positives are ${b.dataset.bbb_pos} experimentally validated peptides.
        Every metric here is a small-sample estimate; the similarity split makes it honest,
        not precise.</li>
        <li>Negatives are random and CPP peptides, so the model separates BBB-penetrating
        peptides from <em>those</em> decoys — not from a natural, drug-like non-BBB library.</li>
        <li>Cell penetration is trained on CPP-vs-random data as a membrane-permeability
        proxy; it says nothing about tumour selectivity.</li>
      </ul>
    </section>`;

  groupedBars(host.querySelector("#bars-a")!, {
    groups: groupsA,
    seriesNames: ["B3Pred published RF", "ours (their split)", "ours (novel-only)"],
    max: 1,
    yLabel: "AUROC",
    annotations,
  });
  groupedBars(host.querySelector("#bars-b")!, {
    groups: groupsB,
    seriesNames: ["AUROC", "AUPRC"],
    max: 1,
    yLabel: "score",
  });
  lineChart(host.querySelector("#line-abstain")!, {
    series: [
      { name: "BBB · similarity split", points: sweep("bbb", "similarity") },
      { name: "BBB · random split", points: sweep("bbb", "random"), dash: true },
      { name: "Cell · similarity split", points: sweep("cpp", "similarity") },
    ],
    xLabel: "abstention rate",
    yLabel: "error among answered",
    xDomain: [0, 0.7],
    yDomain: [0, 0.3],
  });
  lineChart(host.querySelector("#line-reliability")!, {
    series: [
      {
        name: "BBB · similarity",
        points: b.tasks.bbb.similarity.reliability.map((r) => ({ x: r.predicted, y: r.observed })),
      },
      {
        name: "Cell · similarity",
        points: b.tasks.cpp.similarity.reliability.map((r) => ({ x: r.predicted, y: r.observed })),
        dash: true,
      },
    ],
    xLabel: "predicted probability",
    yLabel: "observed frequency",
    xDomain: [0, 1],
    yDomain: [0, 1],
    diagonal: true,
  });
}
