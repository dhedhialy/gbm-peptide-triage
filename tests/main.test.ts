// @vitest-environment jsdom
// Boots src/main.ts the way the browser does: stub fetch, mount #root, then
// drive the tabs and the dial. This is the only test that covers wiring.

import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import modelRaw from "../public/model.json";
import benchRaw from "../public/benchmark.json";

type FetchLike = (input: unknown) => Promise<{ ok: boolean; status: number; json: () => Promise<unknown> }>;

function stubFetch(overrides: Record<string, unknown> = {}): { fail: string | null } {
  const state = { fail: null as string | null };
  const table: Record<string, unknown> = {
    "model.json": modelRaw,
    "benchmark.json": benchRaw,
    ...overrides,
  };
  vi.stubGlobal("fetch", ((input: unknown) => {
    const url = String(input);
    const key = Object.keys(table).find((k) => url.endsWith(k));
    if (state.fail && url.endsWith(state.fail)) {
      return Promise.resolve({ ok: false, status: 404, json: async () => ({}) });
    }
    if (!key) return Promise.reject(new Error(`unexpected fetch ${url}`));
    return Promise.resolve({ ok: true, status: 200, json: async () => table[key] });
  }) as FetchLike);
  return state;
}

async function boot(): Promise<HTMLElement> {
  document.body.innerHTML = `<div id="root"></div>`;
  vi.resetModules();
  await import("../src/main");
  for (let i = 0; i < 4; i++) await new Promise((r) => setTimeout(r, 0));
  return document.getElementById("root")!;
}

beforeEach(() => {
  vi.unstubAllGlobals();
});

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("app bootstrap", () => {
  it("loads the model, auto-analyses the example batch and renders the predict tab", async () => {
    stubFetch();
    const root = await boot();
    expect(root.querySelectorAll("[data-tab]")).toHaveLength(4);
    expect(root.querySelector("[data-tab].active")!.textContent).toContain("Predict");
    expect(root.querySelector<HTMLTextAreaElement>("#batch")!.value).toContain("Angiopep-2");
    expect(root.querySelectorAll("table.results tbody tr.row")).toHaveLength(9);
    expect(root.querySelector<HTMLInputElement>("#dial")!.value).toBe("4");
    expect(root.querySelector("#dial-coverage")!.textContent).toBe("90%");
  });

  it("re-runs conformal sets when the dial moves and keeps the dial note honest", async () => {
    stubFetch();
    const root = await boot();
    const note = root.querySelector("#dial-note")!.textContent!;
    expect(note).toMatch(/measured on held-out clusters at α=0\.1:/);
    expect(note).toMatch(/leaves 39\% unanswered/);
    expect(note).toMatch(/Across the whole dial/);
    const abstainCount = () =>
      Array.from(root.querySelectorAll(".chip")).filter((c) => c.textContent === "abstain").length;

    const dial = root.querySelector<HTMLInputElement>("#dial")!;
    dial.value = "0"; // alpha 0.40, looser coverage target
    dial.dispatchEvent(new Event("input", { bubbles: true }));
    expect(root.querySelector("#dial-coverage")!.textContent).toBe("60%");
    const loosened = abstainCount();
    expect(root.querySelector("#dial-note")!.textContent).toMatch(/at α=0\.4: BBB leaves 32% unanswered/);

    dial.value = "5"; // alpha 0.05, tighter coverage target
    dial.dispatchEvent(new Event("input", { bubbles: true }));
    expect(root.querySelector("#dial-coverage")!.textContent).toBe("95%");
    expect(root.querySelector("#dial-note")!.textContent).toMatch(/at α=0\.05: BBB leaves 60% unanswered/);
    expect(abstainCount()).toBeGreaterThanOrEqual(loosened);
    expect(root.querySelectorAll("table.results tbody tr.row")).toHaveLength(9);
  });

  it("switches tabs and keeps the batch across them", async () => {
    stubFetch();
    const root = await boot();
    const clickTab = (id: string) =>
      root.querySelector<HTMLButtonElement>(`[data-tab="${id}"]`)!.click();

    clickTab("pareto");
    expect(root.querySelector("#panel svg.chart")).not.toBeNull();
    expect(root.querySelectorAll("g.pt").length).toBe(9);
    expect(root.querySelector("[data-tab].active")!.textContent).toContain("Pareto");

    clickTab("validation");
    expect(root.querySelectorAll("#panel svg.chart")).toHaveLength(4);

    clickTab("methods");
    expect(root.querySelector("#panel")!.textContent).toContain("Known failure cases");

    clickTab("predict");
    expect(root.querySelector<HTMLTextAreaElement>("#batch")!.value).toContain("Angiopep-2");
  });

  it("re-analyses an edited batch from the textarea", async () => {
    stubFetch();
    const root = await boot();
    const ta = root.querySelector<HTMLTextAreaElement>("#batch")!;
    ta.value = ">only\nGRKKRRQRRRPPQ";
    ta.dispatchEvent(new Event("input", { bubbles: true }));
    root.querySelector<HTMLButtonElement>("#analyze")!.click();
    expect(root.querySelectorAll("table.results tbody tr.row")).toHaveLength(1);
    expect(root.querySelector(".seq")!.textContent).toBe("GRKKRRQRRRPPQ");
  });

  it("explains a missing model instead of rendering a broken page", async () => {
    stubFetch().fail = "model.json";
    const root = await boot();
    expect(root.textContent).toContain("Could not load the model");
    expect(root.textContent).toContain("scripts/train.py");
  });
});
