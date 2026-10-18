// @vitest-environment jsdom
// The UI is the deliverable, so it is exercised here the way a user drives it:
// render a tab into a real DOM, click the controls, assert the output.

import { beforeEach, describe, expect, it, vi } from "vitest";
import modelRaw from "../public/model.json";
import benchRaw from "../public/benchmark.json";
import type { ModelFile } from "../src/lib/model";
import { EXAMPLE_BATCH, analyzeBatch } from "../src/lib/triage";
import { parseInput } from "../src/lib/parse";
import type { Actions, AppState, Benchmark, TabId } from "../src/ui/context";
import { renderPredict } from "../src/ui/predict";
import { renderPareto } from "../src/ui/pareto";
import { renderValidation } from "../src/ui/validation";
import { renderMethods } from "../src/ui/methods";

const model = modelRaw as unknown as ModelFile;
const benchmark = benchRaw as unknown as Benchmark;

function makeState(overrides: Partial<AppState> = {}): AppState {
  return {
    model,
    benchmark,
    input: EXAMPLE_BATCH,
    alphaIndex: 4,
    results: [],
    selected: null,
    picked: null,
    paretoX: "bbb",
    paretoY: "cpp",
    showGhosts: true,
    ...overrides,
  };
}

function makeHost(): HTMLElement {
  const host = document.createElement("div");
  document.body.appendChild(host);
  return host;
}

function makeActions() {
  const calls: Record<string, unknown[]> = {
    analyze: [], select: [], pick: [], axis: [], ghosts: [], goTab: [],
  };
  const actions: Actions = {
    setInput: vi.fn(),
    analyze: vi.fn(() => calls.analyze.push(true)),
    setAlpha: vi.fn(),
    select: vi.fn((i) => calls.select.push(i)),
    pick: vi.fn((i) => calls.pick.push(i)),
    setAxis: vi.fn((a, k) => calls.axis.push([a, k])),
    toggleGhosts: vi.fn(() => calls.ghosts.push(true)),
    reRender: vi.fn(),
    goTab: vi.fn((t: TabId) => calls.goTab.push(t)),
  };
  return { actions, calls };
}

beforeEach(() => {
  document.body.innerHTML = "";
});

describe("predict tab", () => {
  it("renders an empty state before a batch is analyzed", () => {
    const host = makeHost();
    const { actions } = makeActions();
    renderPredict(host, makeState(), actions);
    expect(host.querySelector<HTMLTextAreaElement>("#batch")!.value).toBe(EXAMPLE_BATCH);
    expect(host.querySelectorAll("table.results tbody tr")).toHaveLength(0);
    expect(host.querySelector(".empty")).not.toBeNull();
  });

  it("fills the triage table from results and wires the analyze/example buttons", () => {
    const host = makeHost();
    const { actions, calls } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPredict(host, makeState({ results }), actions);

    expect(host.querySelectorAll("table.results tbody tr.row")).toHaveLength(9);
    expect(host.querySelectorAll(".chip-pos, .chip-neg, .chip-abstain, .chip-ood").length).toBeGreaterThan(0);
    expect(host.querySelectorAll(".bar-mini")).toHaveLength(9);
    expect(host.querySelector(".stats")!.textContent).toContain("9");

    host.querySelector<HTMLButtonElement>("#analyze")!.click();
    expect(calls.analyze).toHaveLength(1);
    host.querySelector<HTMLButtonElement>("#example")!.click();
    expect(actions.setInput).toHaveBeenCalledWith(EXAMPLE_BATCH);
    expect(calls.analyze).toHaveLength(2);
  });

  it("shows reasons, attribution, neighbours and the heuristic breakdown when a row is opened", () => {
    const host = makeHost();
    const { actions } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPredict(host, makeState({ results, selected: 0 }), actions);

    const detail = host.querySelector(".detail-grid")!;
    expect(detail).not.toBeNull();
    expect(detail.querySelectorAll(".reasons li").length).toBeGreaterThanOrEqual(3);
    expect(detail.querySelectorAll(".attr-row").length).toBeGreaterThan(0);
    expect(detail.querySelectorAll("table.mini-table").length).toBeGreaterThanOrEqual(3);
    expect(detail.textContent).toContain("heuristic index");
    expect(detail.textContent).toContain("Forest attribution");
  });

  it("routes row toggles and invalid peptides", () => {
    const host = makeHost();
    const { actions, calls } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPredict(host, makeState({ results }), actions);
    host.querySelector<HTMLButtonElement>("[data-open]")!.click();
    expect(calls.select).toEqual([0]);

    const bad = analyzeBatch(model, parseInput("PEPTIDE1"), 4);
    const host2 = makeHost();
    renderPredict(host2, makeState({ input: "PEPTIDE1", results: bad }), actions);
    expect(host2.querySelector(".row-error")!.textContent).toContain("non-standard");
  });
});

describe("pareto tab", () => {
  it("plots the batch, the reference cloud and the front, and reports picks", () => {
    const host = makeHost();
    const { actions, calls } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPareto(host, makeState({ results }), actions);

    const svg = host.querySelector("svg.chart")!;
    expect(svg).not.toBeNull();
    expect(svg.querySelectorAll("circle.pt-ghost").length).toBeGreaterThan(20);
    expect(host.querySelectorAll("g.pt").length).toBe(9);
    expect(host.querySelectorAll("g.pt-front").length).toBeGreaterThan(0);
    expect(host.querySelector(".chips")!.querySelectorAll("button.chip-pick").length).toBeGreaterThan(0);

    const pt = host.querySelector<SVGGElement>("g.pt")!;
    pt.dispatchEvent(new MouseEvent("click", { bubbles: true }));
    expect(calls.pick).toEqual([0]);

    host.querySelector<HTMLButtonElement>("[data-pick]")!.click();
    expect(calls.pick.length).toBe(2);
  });

  it("re-plots when an axis changes and sends the axis back to the app", () => {
    const host = makeHost();
    const { actions, calls } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPareto(host, makeState({ results }), actions);
    const sel = host.querySelector<HTMLSelectElement>("#px")!;
    sel.value = "immune";
    sel.dispatchEvent(new Event("change", { bubbles: true }));
    expect(calls.axis).toEqual([["x", "immune"]]);
  });

  it("explains itself when there is no batch yet", () => {
    const host = makeHost();
    const { actions, calls } = makeActions();
    renderPareto(host, makeState(), actions);
    expect(host.querySelector(".empty")!.textContent).toContain("Predict");
    host.querySelector<HTMLButtonElement>("[data-goto]")!.click();
    expect(calls.goTab).toEqual(["predict"]);
  });

  it("keeps out-of-domain peptides out of the front but still shows them", () => {
    const host = makeHost();
    const { actions } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    const oodIndex = results.findIndex((r) => r.bbb!.ood || r.cpp!.ood);
    expect(oodIndex).toBeGreaterThanOrEqual(0);
    renderPareto(host, makeState({ results }), actions);
    const frontLabels = Array.from(host.querySelectorAll("button.chip-pick")).map((b) => b.textContent);
    expect(frontLabels).not.toContain(results[oodIndex].label);
    expect(host.querySelectorAll("g.pt")).toHaveLength(9);
    expect(host.textContent).toContain("excluded from ranking");
  });

  it("pins a peptide and names its dominators", () => {
    const host = makeHost();
    const { actions } = makeActions();
    const results = analyzeBatch(model, parseInput(EXAMPLE_BATCH), 4);
    renderPareto(host, makeState({ results, picked: 5 }), actions);
    expect(host.querySelector(".picked")!.textContent).toContain(results[5].label);
  });
});

describe("validation tab", () => {
  it("draws the B3Pred audit, the protocol comparison and the dial curves", () => {
    const host = makeHost();
    renderValidation(host, makeState());
    expect(host.querySelectorAll("svg.chart")).toHaveLength(4);
    expect(host.querySelectorAll(".legend").length).toBe(4);
    expect(host.querySelectorAll("table.mini-table.wide").length).toBe(2);
    expect(host.textContent).toContain("B3Pred");
    expect(host.textContent).toContain((benchmark.tasks.bbb.similarity.n_clusters ?? 0).toString());
    expect(host.textContent).toContain(benchmark.b3pred_published_splits.datasets.d3.our_auprc.toFixed(2));
  });

  it("shows a loading note until the benchmark arrives", () => {
    const host = makeHost();
    renderValidation(host, makeState({ benchmark: null }));
    expect(host.querySelector(".empty")!.textContent).toContain("Loading");
  });
});

describe("methods tab", () => {
  it("lists the descriptor set, the conformal recipe and the limitations", () => {
    const host = makeHost();
    renderMethods(host, model);
    expect(host.textContent).toContain("frac_V");
    expect(host.textContent).toContain("34 transparent descriptors");
    expect(host.textContent).toContain("LAC nonconformity");
    expect(host.textContent).toContain("Known failure cases");
    expect(host.textContent).toContain("Pharmaceutics 2021");
  });
});
