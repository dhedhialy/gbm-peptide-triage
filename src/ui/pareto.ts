import { scatter, type ScatterPoint } from "../viz/charts";
import { paretoFront } from "../lib/pareto";
import { OBJECTIVES, analyzeBatch, type ObjectiveKey, type PeptideResult } from "../lib/triage";
import type { AppState, Actions } from "./context";
import { esc } from "./dom";

// reference cloud: a sample of training positives, scored once and cached
let ghostObjectives: Record<ObjectiveKey, number>[] | null = null;

function ghosts(state: AppState): Record<ObjectiveKey, number>[] {
  if (ghostObjectives) return ghostObjectives;
  const positives = state.model.tasks.bbb.train_seqs.filter((t) => t.y === 1);
  const step = Math.max(1, Math.floor(positives.length / 120));
  const sample = positives.filter((_, i) => i % step === 0);
  const results = analyzeBatch(
    state.model,
    sample.map((t, i) => ({ label: `ref ${i}`, seq: t.s })),
    4, // α = 0.10
  );
  ghostObjectives = results
    .filter((r) => r.objectives)
    .map((r) => r.objectives!);
  return ghostObjectives;
}

export function renderPareto(host: HTMLElement, state: AppState, actions: Actions): void {
  const valid: PeptideResult[] = state.results.filter((r) => !r.error && r.objectives);
  // Out-of-domain peptides keep their heuristic scores (rules, not learned
  // weights) but their learned objectives are not trustworthy, so they are
  // plotted as unranked and never enter the front.
  const rankable = valid.map((r) => !r.bbb!.ood && !r.cpp!.ood);
  const label = (k: ObjectiveKey) => OBJECTIVES.find((o) => o.key === k)!.label;
  const axisOptions = (sel: ObjectiveKey) =>
    OBJECTIVES.map(
      (o) => `<option value="${o.key}"${o.key === sel ? " selected" : ""}>${o.label} <small>(${o.source})</small></option>`,
    ).join("");

  const controls = `
    <div class="row row-between wrap">
      <div class="row wrap">
        <label class="inline">x-axis <select id="px">${axisOptions(state.paretoX)}</select></label>
        <label class="inline">y-axis <select id="py">${axisOptions(state.paretoY)}</select></label>
        <label class="inline check"><input type="checkbox" id="ghosts"${state.showGhosts ? " checked" : ""}>
          reference cloud <span class="tag">training BBB+</span></label>
      </div>
      <span class="hint">hover to inspect · click a point to pin</span>
    </div>`;

  if (valid.length === 0) {
    host.innerHTML = `${controls}<div class="empty">
      <p>The Pareto view plots the batch from the <strong>Predict</strong> tab.
      Analyze at least two peptides, then come back — or
      <button class="link" data-goto="predict">go to Predict</button>.</p>
    </div>`;
    wireControls(host, actions);
    return;
  }

  const dimsOf = (r: PeptideResult) => [
    r.objectives!.bbb, r.objectives!.cpp, r.objectives!.immune, r.objectives!.mfg,
  ];
  const rankableIndexes = rankable.map((ok, i) => (ok ? i : -1)).filter((i) => i >= 0);
  const front2full = paretoFront(
    rankableIndexes.map((i) => ({
      dims: [valid[i].objectives![state.paretoX], valid[i].objectives![state.paretoY]],
    })),
  );
  const front2 = valid.map(() => false);
  rankableIndexes.forEach((i, k) => {
    front2[i] = front2full[k];
  });
  const front4full = paretoFront(rankableIndexes.map((i) => ({ dims: dimsOf(valid[i]) })));
  const front4 = valid.map(() => false);
  rankableIndexes.forEach((i, k) => {
    front4[i] = front4full[k];
  });

  const pts: ScatterPoint[] = [];
  if (state.showGhosts) {
    for (const g of ghosts(state)) {
      pts.push({ x: g[state.paretoX], y: g[state.paretoY], label: "reference BBB+", sub: "", kind: "ghost" });
    }
  }
  const start = pts.length;
  valid.forEach((r, i) => {
    const o = r.objectives!;
    const onFront = front2[i];
    pts.push({
      x: o[state.paretoX],
      y: o[state.paretoY],
      label: esc(r.label),
      sub:
        `${esc(r.seq)}<br>${OBJECTIVES.map((ob) => `${ob.short} ${o[ob.key].toFixed(0)}`).join(" · ")}` +
        `${front4[i] ? "<br><em>on the 4-objective front</em>" : ""}` +
        `${!rankable[i] ? "<br><em>out-of-domain: not ranked</em>" : ""}`,
      kind: onFront ? "front" : !rankable[i] || r.bbb!.verdict === "abstain" ? "abstain" : "plain",
    });
  });

  const frontList = valid
    .map((r, i) => (front2[i] ? { r, i } : null))
    .filter(Boolean) as { r: PeptideResult; i: number }[];
  const picked = state.picked !== null && state.picked < valid.length ? state.picked : null;
  const dominators =
    picked === null || !rankable[picked]
      ? []
      : valid
          .map((r, i) => ({ r, i }))
          .filter(
            ({ i }) =>
              i !== picked && rankable[i] &&
              dominates(
                [valid[i].objectives![state.paretoX], valid[i].objectives![state.paretoY]],
                [valid[picked].objectives![state.paretoX], valid[picked].objectives![state.paretoY]],
              ),
          )
          .slice(0, 4);

  host.innerHTML = `
    ${controls}
    <div class="grid-pareto">
      <div class="card" id="scatter-host"></div>
      <div class="card">
        <h3>Front on the selected axes
          <span class="tag">${frontList.length} of ${valid.length}</span></h3>
        <p class="hint">A peptide is on the front when no other peptide is at least as good on
        both axes and strictly better on one. All four objectives run in the background:
        <strong>${front4.filter(Boolean).length} of ${valid.length}</strong> peptides are
        non-dominated under all four.${
          rankable.includes(false)
            ? ` Out-of-domain peptides (grey, dashed) are plotted but excluded from ranking — ` +
              `their learned objectives are not trustworthy; their heuristic scores still are.`
            : ""
        }</p>
        <div class="chips">
          ${frontList
            .map(
              ({ r, i }) =>
                `<button class="chip chip-pick${i === picked ? " active" : ""}" data-pick="${i}">${esc(r.label)}</button>`,
            )
            .join("") || `<span class="hint">no front members</span>`}
        </div>
        ${
          picked === null
            ? `<p class="hint">Pin a point to see its trade-offs here.</p>`
            : `<div class="picked">
                <h3>${esc(valid[picked].label)}</h3>
                <p class="mono">${esc(valid[picked].seq)}</p>
                <table class="mini-table">${OBJECTIVES.map(
                  (o) =>
                    `<tr><td>${o.label}</td><td class="num">${valid[picked].objectives![o.key].toFixed(0)}</td>
                     <td class="hint">${o.source}</td></tr>`,
                ).join("")}</table>
                <p class="hint">${
                  !rankable[picked]
                    ? "Out-of-domain: shown for inspection, excluded from the front."
                    : dominators.length === 0
                    ? "Nothing in the batch beats it on both selected axes."
                    : `On these axes it is dominated by ${dominators
                        .map(({ r }) => `<strong>${esc(r.label)}</strong>`)
                        .join(", ")}.`
                }</p>
              </div>`
        }
      </div>
    </div>`;

  scatter(host.querySelector("#scatter-host")!, {
    points: pts,
    xLabel: label(state.paretoX),
    yLabel: label(state.paretoY),
    xDomain: [0, 100],
    yDomain: [0, 100],
    picked: picked === null ? undefined : start + picked,
    onPick: (i) => {
      if (i < start) return;
      const idx = i - start;
      actions.pick(state.picked === idx ? null : idx);
    },
  });

  wireControls(host, actions);
  host.querySelectorAll<HTMLElement>("[data-pick]").forEach((el) =>
    el.addEventListener("click", () => actions.pick(Number(el.dataset.pick))),
  );
}

function wireControls(host: HTMLElement, actions: Actions): void {
  host.querySelector<HTMLSelectElement>("#px")?.addEventListener("change", (e) =>
    actions.setAxis("x", (e.target as HTMLSelectElement).value as ObjectiveKey),
  );
  host.querySelector<HTMLSelectElement>("#py")?.addEventListener("change", (e) =>
    actions.setAxis("y", (e.target as HTMLSelectElement).value as ObjectiveKey),
  );
  host.querySelector<HTMLInputElement>("#ghosts")?.addEventListener("change", () =>
    actions.toggleGhosts(),
  );
  host.querySelectorAll<HTMLElement>("[data-goto]").forEach((el) =>
    el.addEventListener("click", () => actions.goTab("predict")),
  );
}

function dominates(a: number[], b: number[]): boolean {
  let strictly = false;
  for (let k = 0; k < a.length; k++) {
    if (a[k] < b[k] - 1e-9) return false;
    if (a[k] > b[k] + 1e-9) strictly = true;
  }
  return strictly;
}
