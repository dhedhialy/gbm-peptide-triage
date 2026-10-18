import type { ModelFile } from "../lib/model";
import type { PeptideResult, ObjectiveKey } from "../lib/triage";

export interface ProtocolStats {
  auroc_mean: number;
  auroc_sd: number;
  auprc_mean: number;
  auprc_sd: number;
  n_seeds: number;
  q: number;
  n_train: number;
  n_cal: number;
  n_test: number;
  n_clusters?: number;
  conformal: {
    coverage: number;
    abstain_rate: number;
    empty_rate: number;
    singleton_acc: number;
    mean_set_size: number;
  };
  reliability: { predicted: number; observed: number; count: number }[];
  alpha_sweep: { alpha: number; abstain_rate: number; error_rate: number; coverage: number }[];
  test_to_train_sim: { p10: number; median: number; p90: number; frac_below_ood: number };
}

export interface Benchmark {
  alpha: number;
  model: Record<string, unknown>;
  similarity: { definition: string; cluster_threshold: number; ood_threshold: number };
  dataset: { bbb_pos: number; cpp_pos: number; random_neg: number };
  tasks: Record<"bbb" | "cpp", { random: ProtocolStats; similarity: ProtocolStats }>;
  b3pred_published_splits: {
    citation: string;
    note: string;
    datasets: Record<
      string,
      {
        published_rf_val_auroc: number;
        our_auroc: number;
        our_auprc: number;
        n_train: number;
        n_valid: number;
        "frac_valid_ge_0.8_similar_to_train": number;
        "frac_valid_ge_0.6_similar_to_train": number;
        median_max_sim: number;
        our_auroc_novel_only: number | null;
        n_novel_only: number;
      }
    >;
  };
}

export type TabId = "predict" | "pareto" | "validation" | "methods";

export interface AppState {
  model: ModelFile;
  benchmark: Benchmark | null;
  input: string;
  alphaIndex: number;
  results: PeptideResult[];
  selected: number | null;
  picked: number | null;
  paretoX: ObjectiveKey;
  paretoY: ObjectiveKey;
  showGhosts: boolean;
}

export interface Actions {
  setInput(text: string): void;
  analyze(): void;
  setAlpha(i: number): void;
  select(i: number | null): void;
  pick(i: number | null): void;
  setAxis(axis: "x" | "y", key: ObjectiveKey): void;
  toggleGhosts(): void;
  reRender(): void;
  goTab(tab: TabId): void;
}
