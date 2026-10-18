// Batch triage: runs both conformal predictors plus the two heuristic indexes
// and assembles the four objectives used by the Pareto view.

import { analyzeTask, type ModelFile, type TaskResult } from "./model";
import { immuneCompatibility, manufacturability, type HeuristicResult } from "./heuristics";
import type { ParsedPeptide } from "./parse";

export const ALPHAS = [0.40, 0.30, 0.20, 0.15, 0.10, 0.05];

export type ObjectiveKey = "bbb" | "cpp" | "immune" | "mfg";

export const OBJECTIVES: { key: ObjectiveKey; label: string; short: string; source: string }[] = [
  { key: "bbb", label: "BBB penetration", short: "BBB", source: "trained, conformal" },
  { key: "cpp", label: "Cell penetration", short: "Cell", source: "trained, conformal" },
  { key: "immune", label: "Immune compatibility", short: "Immune", source: "heuristic index" },
  { key: "mfg", label: "Manufacturability", short: "Mfg", source: "heuristic index" },
];

export interface PeptideResult {
  label: string;
  seq: string;
  error?: string;
  bbb?: TaskResult;
  cpp?: TaskResult;
  immune?: HeuristicResult;
  mfg?: HeuristicResult;
  objectives?: Record<ObjectiveKey, number>;
}

export function analyzeBatch(
  model: ModelFile,
  peptides: ParsedPeptide[],
  alphaIndex: number,
): PeptideResult[] {
  const alpha = ALPHAS[alphaIndex];
  const key = String(alpha);
  const coverage = 1 - alpha;
  return peptides.map((p) => {
    if (p.error) return { label: p.label, seq: p.seq, error: p.error };
    const bbb = analyzeTask(
      model.tasks.bbb, p.seq, model.tasks.bbb.qs[key] ?? model.tasks.bbb.q,
      model.ood_sim, model.features, "BBB-penetrating", coverage,
    );
    const cpp = analyzeTask(
      model.tasks.cpp, p.seq, model.tasks.cpp.qs[key] ?? model.tasks.cpp.q,
      model.ood_sim, model.features, "cell-penetrating", coverage,
    );
    const immune = immuneCompatibility(p.seq);
    const mfg = manufacturability(p.seq);
    return {
      label: p.label,
      seq: p.seq,
      bbb,
      cpp,
      immune,
      mfg,
      objectives: {
        bbb: bbb.prob * 100,
        cpp: cpp.prob * 100,
        immune: immune.score,
        mfg: mfg.score,
      },
    };
  });
}

export const EXAMPLE_BATCH = `>Angiopep-2 (literature BBB shuttle)
TFFYGGSRGKRNNFKLKFD
>GE11 (EGFR-targeting, not a BBB shuttle)
YHWGYTPQNVI
>TAT (cell-penetrating, BBB-poor)
GRKKRRQRRRPPQ
>Penetratin (cell-penetrating, BBB-poor)
RQIKIWFQNRRMKWKK
>TGN peptide (B3Pdb validated BBB+)
TGNYKALHPHNG
>B3Pdb positive
CDIFTNSRGKRA
>Scrambled Angiopep-2 (designed control)
DKFKLNNRKGGSGYFFT
>B3Pred CPP negative
LLRARWRRRRSRRFR
>33-mer (outside the 5–30 aa training range: watch it get flagged out-of-domain)
GAVLIPKGAAQTGAVLIPKGAAQTGAVLIPKGA`;
