// Forest inference, conformal prediction sets, and local attribution.
// The forest encoding is produced by export_forest() in scripts/train.py.

import { features } from "./features";
import { nearestNeighbors, maxSimTo, type Neighbor } from "./identity";

export interface Tree {
  left: number[];
  right: number[];
  feature: number[];
  thr: number[];
  prob: number[];
}

export interface TaskModel {
  n_trees: number;
  trees: Tree[];
  importance: number[];
  q: number;
  qs: Record<string, number>;
  train_seqs: { s: string; y: number }[];
  train_len_range: [number, number];
  n_cal: number;
  n_test: number;
  test_indices: number[];
}

export interface Fixture {
  task: string;
  s: string;
  features: number[];
  prob: number;
  set: [boolean, boolean];
  sim_to_train: number;
  top_attr: [number, number][];
}

export interface ModelFile {
  version: number;
  features: string[];
  alpha: number;
  ood_sim: number;
  tasks: { bbb: TaskModel; cpp: TaskModel };
  fixtures: Fixture[];
}

export function forestProb(task: TaskModel, x: number[]): number {
  let total = 0;
  for (const t of task.trees) {
    let i = 0;
    while (t.left[i] !== -1) {
      i = x[t.feature[i]] <= t.thr[i] ? t.left[i] : t.right[i];
    }
    total += t.prob[i];
  }
  return total / task.n_trees;
}

/** Per-feature probability deltas along each tree path, averaged (root-to-leaf). */
export function attribution(task: TaskModel, x: number[]): [number, number][] {
  const acc = new Map<number, number>();
  for (const t of task.trees) {
    let i = 0;
    while (t.left[i] !== -1) {
      const j = x[t.feature[i]] <= t.thr[i] ? t.left[i] : t.right[i];
      const f = t.feature[i];
      acc.set(f, (acc.get(f) ?? 0) + (t.prob[j] - t.prob[i]));
      i = j;
    }
  }
  return [...acc.entries()]
    .map(([f, d]) => [f, d / task.n_trees] as [number, number])
    .filter(([, d]) => Math.abs(d) > 1e-9)
    .sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
}

/** LAC split-conformal set: include label y when 1 - p(y) <= q. */
export function conformalSet(prob: number, q: number): [boolean, boolean] {
  const eps = 1e-12;
  return [prob <= q + eps, 1 - prob <= q + eps];
}

export type Verdict = "positive" | "negative" | "abstain";

export interface TaskResult {
  prob: number;
  set: [boolean, boolean];
  q: number;
  verdict: Verdict;
  ood: boolean;
  maxSim: number;
  neighbors: Neighbor[];
  topAttr: { feature: string; delta: number }[];
  reasons: string[];
}

export function verdictOf(set: [boolean, boolean]): Verdict {
  if (set[0] && set[1]) return "abstain";
  if (set[0]) return "negative";
  if (set[1]) return "positive";
  return "abstain"; // empty set: no label survives the test
}

export function analyzeTask(
  task: TaskModel,
  seq: string,
  q: number,
  oodSimThreshold: number,
  featureNames: string[],
  posLabel: string,
  coverage: number,
): TaskResult {
  const x = features(seq);
  const prob = forestProb(task, x);
  const set = conformalSet(prob, q);
  const nn = nearestNeighbors(seq, task.train_seqs, 3);
  const maxSim = maxSimTo(seq, task.train_seqs)?.sim ?? 0;
  const [lo, hi] = task.train_len_range;
  const lengthOut = seq.length < lo || seq.length > hi;
  const simLow = maxSim < oodSimThreshold;
  const ood = simLow || lengthOut;
  const verdict = verdictOf(set);
  const topAttr = attribution(task, x)
    .slice(0, 4)
    .map(([f, d]) => ({ feature: featureNames[f], delta: d }));
  const cov = `${(100 * coverage).toFixed(0)}%`;
  const reasons: string[] = [];
  if (simLow) {
    reasons.push(
      `Outside the training distribution: the closest training peptide is only ` +
        `${(maxSim * 100).toFixed(0)}% similar (flag below the ${(oodSimThreshold * 100).toFixed(0)}% ` +
        `gate, which sits at the bottom decile of held-out peptide similarity). ` +
        `The ${cov} coverage guarantee does not hold here.`,
    );
  }
  if (lengthOut) {
    reasons.push(
      `Length ${seq.length} is outside the ${lo}–${hi} aa range seen in training, so every ` +
        `descriptor extrapolates. Treated as out-of-domain.`,
    );
  }
  if (ood) {
    reasons.push(
      `Verdict withheld: this is an out-of-domain flag, and a conformal guarantee is void ` +
        `outside the distribution it was calibrated on. The probability below is the raw ` +
        `ensemble output, shown for inspection only.`,
    );
  } else if (verdict === "abstain") {
    reasons.push(
      `Both labels survive the conformal test at the shipped ${cov} coverage target: ` +
        `the model cannot separate them for this sequence.`,
    );
  }
  if (!ood && verdict !== "abstain") {
    reasons.push(
      `Conformal singleton {${verdict === "positive" ? posLabel : "not " + posLabel}}; ` +
        `${cov} marginal coverage on in-distribution peptides.`,
    );
  }
  if (nn.length) {
    reasons.push(
      `Closest known peptides: ${nn
        .map((n) => `${n.seq} (${n.label ? posLabel : "negative"}, ${(n.sim * 100).toFixed(0)}%)`)
        .join("; ")}.`,
    );
  }
  return { prob, set, q, verdict, ood, maxSim, neighbors: nn, topAttr, reasons };
}
