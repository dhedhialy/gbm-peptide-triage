// Transparent heuristic indexes. These are NOT trained or calibrated models —
// they are deterministic indexes over published scales and liability rules,
// and the UI labels them as such. Every component is shown with its points.

import { features } from "./features";

export interface Component {
  label: string;
  points: number;
  detail: string;
}

export interface HeuristicResult {
  score: number; // 0-100, higher = better (immune compatibility / manufacturability)
  components: Component[];
}

const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));

// MHC-II binding groove prefers hydrophobic/aromatic anchors at P1 and
// secondary anchors P4/P6/P9 (open-ended groove, 9-mer core).
const ANCHOR = new Set("AVILMFWY");

function kmers(seq: string, k: number): string[] {
  const out: string[] = [];
  for (let i = 0; i + k <= seq.length; i++) out.push(seq.slice(i, i + k));
  return out;
}

/** T-cell presentation propensity (higher score = lower risk). */
export function immuneCompatibility(seq: string): HeuristicResult {
  const components: Component[] = [];
  const frames = kmers(seq, 9);
  let meanAnchor: number;
  if (frames.length > 0) {
    let hits = 0;
    for (const f of frames) {
      const pos = [0, 3, 5, 8];
      hits += pos.filter((p) => ANCHOR.has(f[p])).length / pos.length;
    }
    meanAnchor = hits / frames.length;
  } else {
    // shorter than one 9-mer: weaker class-II competence, partial anchor read
    const a = [...seq].filter((c) => ANCHOR.has(c)).length / seq.length;
    meanAnchor = a * 0.5;
  }
  const anchorPoints = 75 * meanAnchor;
  components.push({
    label: "MHC-II anchor density",
    points: -anchorPoints, // risk
    detail: `${(meanAnchor * 100).toFixed(0)}% of 9-mer frames carry a hydrophobic P1/P4/P6/P9 anchor pattern`,
  });
  const lenFactor =
    seq.length >= 13 && seq.length <= 25 ? 1 : seq.length >= 9 ? 0.5 : 0.15;
  const lenPoints = 25 * lenFactor;
  components.push({
    label: "Class-II length window",
    points: -lenPoints, // risk
    detail:
      lenFactor === 1
        ? `${seq.length} aa sits in the 13–25 aa optimal MHC-II range`
        : `${seq.length} aa is outside the optimal 13–25 aa range (factor ${lenFactor})`,
  });
  const score = clamp(100 - (anchorPoints + lenPoints), 0, 100);
  return { score: Math.round(score), components };
}

const DEAMIDATION = /N[GSTHA]/g;
const ISOMERIZATION = /D[GST]/g;

/** Manufacturing/developability index from charge, aggregation and liabilities. */
export function manufacturability(seq: string): HeuristicResult {
  const components: Component[] = [];
  const charge = Math.abs(features(seq)[1]);
  const chargeScore =
    charge >= 2 && charge <= 6
      ? 1
      : charge < 2
        ? charge / 2
        : Math.max(0, 1 - (charge - 6) / 4);
  let run = 0;
  let best = 0;
  for (const c of seq) {
    run = "ACILFVM".includes(c) ? run + 1 : 0;
    best = Math.max(best, run);
  }
  const runScore = best <= 4 ? 1 : best >= 8 ? 0 : 1 - (best - 4) / 4;
  const base = 100 * (0.5 * chargeScore + 0.5 * runScore);
  components.push({
    label: "Solubility (net charge)",
    points: 50 * chargeScore,
    detail: `|net charge| = ${charge.toFixed(1)} at pH 7.4 (2–6 is the favorable band)`,
  });
  components.push({
    label: "Aggregation (hydrophobic run)",
    points: 50 * runScore,
    detail: `longest hydrophobic run = ${best} residue${best === 1 ? "" : "s"}`,
  });

  const penalties: [string, number, number, string][] = [
    ["Deamidation motifs", count(DEAMIDATION, seq), 8, "Asn-Gly/Ser/Thr/His/Ala"],
    ["Oxidation (Met)", count(/M/g, seq), 4, "methionine sulfoxide formation"],
    ["Oxidation (Trp)", count(/W/g, seq), 6, "tryptophan oxidation"],
    ["Free cysteines", count(/C/g, seq), 8, "disulfide scrambling / conjugation"],
    ["Asp isomerization motifs", count(ISOMERIZATION, seq), 6, "Asp-Gly/Ser/Thr"],
  ];
  let totalPenalty = 0;
  for (const [label, n, per, detail] of penalties) {
    if (n === 0) continue;
    const pts = Math.min(n * per, per * 4);
    totalPenalty += pts;
    components.push({
      label,
      points: -pts,
      detail: `${n} occurrence${n === 1 ? "" : "s"} of ${detail}`,
    });
  }
  if (seq[0] === "Q") {
    totalPenalty += 6;
    components.push({ label: "N-terminal Gln", points: -6, detail: "pyroglutamate cyclization risk" });
  }
  const score = clamp(base - totalPenalty, 0, 100);
  return { score: Math.round(score), components };
}

function count(re: RegExp, seq: string): number {
  return (seq.match(re) ?? []).length;
}
