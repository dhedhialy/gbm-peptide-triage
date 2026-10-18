import { EXAMPLE_BATCH, type ObjectiveKey } from "../lib/triage";
import { features } from "../lib/features";
import type { TaskResult } from "../lib/model";
import type { HeuristicResult } from "../lib/heuristics";
import type { AppState, Actions } from "./context";
import { esc, delegate, pct } from "./dom";

let sortKey: "index" | ObjectiveKey = "index";

function chip(t: TaskResult, prefix: string): string {
  if (t.ood) return `<span class="chip chip-ood">OOD</span>`;
  const map = {
    positive: [`chip chip-pos`, `${prefix}+`],
    negative: [`chip chip-neg`, `${prefix}−`],
    abstain: [`chip chip-abstain`, "abstain"],
  } as const;
  const [cls, label] = map[t.verdict];
  return `<span class="${cls}" title="${esc(t.reasons[0] ?? "")}">${label}</span>`;
}

function bar(p: number): string {
  return `<span class="bar-mini"><i style="width:${(p * 100).toFixed(1)}%"></i></span>`;
}

function detailRow(r: AppState["results"][number]): string {
  if (r.error) return "";
  const x = features(r.seq);
  const valueOf = (name: string) => x[state_featureNames.indexOf(name)] ?? 0;
  const attr = r.bbb!.topAttr;
  const maxAbs = Math.max(...attr.map((a) => Math.abs(a.delta)), 1e-6);
  return `
  <tr class="detail-row">
    <td colspan="10">
      <div class="detail-grid">
        <div class="card inner">
          <h3>Why — BBB verdict</h3>
          <p class="prob-line">P(BBB-penetrating) = <strong>${r.bbb!.prob.toFixed(3)}</strong>
            &nbsp;·&nbsp; conformal set <code>{${[
              r.bbb!.set[1] ? "BBB+" : null,
              r.bbb!.set[0] ? "BBB−" : null,
            ].filter(Boolean).join(", ")}}</code>
            &nbsp;·&nbsp; q = ${r.bbb!.q.toFixed(3)} &nbsp;·&nbsp; length ${r.seq.length} aa</p>
          <ul class="reasons">${r.bbb!.reasons.map((s) => `<li>${esc(s)}</li>`).join("")}</ul>
        </div>
        <div class="card inner">
          <h3>Forest attribution (ΔP along the tree paths)</h3>
          <div class="attr">
            ${attr
              .map((a) => {
                const w = (Math.abs(a.delta) / maxAbs) * 50;
                const val = valueOf(a.feature);
                const left = a.delta >= 0 ? 50 : 50 - w;
                return `<div class="attr-row">
                  <span class="attr-name" title="${esc(a.feature)}">${esc(a.feature)}</span>
                  <span class="attr-track"><i class="attr-fill ${a.delta >= 0 ? "pos" : "neg"}" style="left:${left}%;width:${w}%"></i></span>
                  <span class="attr-val">${a.delta >= 0 ? "+" : ""}${(a.delta * 100).toFixed(1)}%</span>
                  <span class="attr-feat">${val.toFixed(2)}</span>
                </div>`;
              })
              .join("")}
          </div>
          <p class="hint">Bars show how much each descriptor moves the probability
          away from the training prior, averaged over the 200 trees. Right of centre = pushes
          toward BBB+.</p>
        </div>
        <div class="card inner">
          <h3>Closest known peptides <span class="tag">training set</span></h3>
          <table class="mini-table">
            <thead><tr><th>sequence</th><th>label</th><th>similarity</th></tr></thead>
            <tbody>${r.bbb!.neighbors
              .map(
                (n) =>
                  `<tr><td class="mono">${esc(n.seq)}</td><td>${n.label ? "BBB+" : "BBB−"}</td><td>${pct(n.sim)}</td></tr>`,
              )
              .join("")}</tbody>
          </table>
          <p class="hint">Similarity = 2·LCS/(len<sub>a</sub>+len<sub>b</sub>).
          Below ${pct(state_oodSim)} the peptide is flagged out-of-domain.</p>
        </div>
        <div class="card inner">
          <h3>Cell-penetration verdict</h3>
          <p class="prob-line">P(cell-penetrating) = <strong>${r.cpp!.prob.toFixed(3)}</strong></p>
          <ul class="reasons">${r.cpp!.reasons.map((s) => `<li>${esc(s)}</li>`).join("")}</ul>
        </div>
        ${heurCard("Immune compatibility", r.immune!, "risk contribution")}
        ${heurCard("Manufacturability", r.mfg!, "penalty")}
      </div>
    </td>
  </tr>`;
}

function heurCard(title: string, h: HeuristicResult, penaltyWord: string): string {
  const rows = h.components
    .filter((c) => Math.abs(c.points) > 0.05)
    .sort((a, b) => a.points - b.points)
    .map(
      (c) => `<tr>
        <td>${esc(c.label)}</td>
        <td class="num ${c.points < 0 ? "neg" : "pos"}">${c.points < 0 ? "−" : "+"}${Math.abs(c.points).toFixed(0)}</td>
        <td class="hint">${esc(c.detail)}${c.points < 0 ? ` (${penaltyWord})` : ""}</td>
      </tr>`,
    )
    .join("");
  return `<div class="card inner">
    <h3>${title} <span class="tag">heuristic index</span></h3>
    <p class="prob-line">score = <strong>${h.score}</strong> / 100</p>
    <table class="mini-table"><tbody>${rows}</tbody></table>
  </div>`;
}

// values injected by render(); module-level so detailRow stays simple
let state_featureNames: string[] = [];
let state_oodSim = 0.4;

export function renderPredict(host: HTMLElement, state: AppState, actions: Actions): void {
  state_featureNames = state.model.features;
  state_oodSim = state.model.ood_sim;
  const valid = state.results.filter((r) => !r.error);
  const abstain = valid.filter((r) => r.bbb!.verdict === "abstain" && !r.bbb!.ood).length;
  const ood = valid.filter((r) => r.bbb!.ood).length;
  const confident = valid.filter((r) => r.bbb!.verdict !== "abstain" && !r.bbb!.ood);

  host.innerHTML = `
    <section class="card">
      <h2>1 · Paste peptides</h2>
      <textarea id="batch" spellcheck="false" placeholder="One peptide per line, or FASTA:\n>name\nTFFYGGSRGKRNNFKLKFD">${esc(state.input)}</textarea>
      <div class="row">
        <button id="analyze" class="btn btn-primary">Analyze <kbd>⌘⏎</kbd></button>
        <button id="example" class="btn">Load examples</button>
        <button id="clear" class="btn">Clear</button>
      </div>
    </section>
    <div class="stats strip">
      ${stat("peptides", String(valid.length))}
      ${stat("confident BBB call", String(confident.length))}
      ${stat("abstained (BBB)", String(abstain))}
      ${stat("out-of-domain", String(ood))}
    </div>

    <section class="card" id="results-card">
      <div class="row row-between">
        <h2>2 · Triage table</h2>
        <label class="inline">sort by
          <select id="sort">
            <option value="index"${sortKey === "index" ? " selected" : ""}>input order</option>
            <option value="bbb"${sortKey === "bbb" ? " selected" : ""}>BBB probability</option>
            <option value="cpp"${sortKey === "cpp" ? " selected" : ""}>cell probability</option>
            <option value="immune"${sortKey === "immune" ? " selected" : ""}>immune compatibility</option>
            <option value="mfg"${sortKey === "mfg" ? " selected" : ""}>manufacturability</option>
          </select>
        </label>
      </div>
      ${state.results.length === 0 ? emptyState() : resultsTable(state)}
    </section>`;

  // ---- events ----
  const ta = host.querySelector<HTMLTextAreaElement>("#batch")!;
  ta.addEventListener("input", () => actions.setInput(ta.value));
  host.querySelector("#analyze")!.addEventListener("click", () => actions.analyze());
  host.querySelector("#example")!.addEventListener("click", () => {
    actions.setInput(EXAMPLE_BATCH);
    actions.analyze();
  });
  host.querySelector("#clear")!.addEventListener("click", () => {
    actions.setInput("");
    actions.analyze();
  });
  host.querySelector<HTMLSelectElement>("#sort")?.addEventListener("change", (e) => {
    sortKey = (e.target as HTMLSelectElement).value as typeof sortKey;
    actions.reRender();
  });
  delegate<HTMLElement>(host, "[data-open]", "click", (el) => {
    const i = Number(el.dataset.open);
    actions.select(state.selected === i ? null : i);
  });
  delegate<HTMLElement>(host, "[data-example]", "click", () => {
    actions.setInput(EXAMPLE_BATCH);
    actions.analyze();
  });
}

function stat(label: string, value: string): string {
  return `<div class="stat"><span class="stat-value">${value}</span><span class="stat-label">${label}</span></div>`;
}

function emptyState(): string {
  return `<div class="empty">
    <p>No batch yet. Paste peptides on the left, or
    <button class="link" data-example>load the example batch</button> —
    literature peptides plus validated B3Pdb positives.</p>
  </div>`;
}

function resultsTable(state: AppState): string {
  const idx = state.results.map((_, i) => i);
  if (sortKey !== "index") {
    idx.sort((a, b) => {
      const av = state.results[a].objectives?.[sortKey as ObjectiveKey] ?? -1;
      const bv = state.results[b].objectives?.[sortKey as ObjectiveKey] ?? -1;
      return bv - av || a - b;
    });
  }
  const rows = idx
    .map((i) => {
      const r = state.results[i];
      if (r.error) {
        return `<tr class="row row-error">
          <td class="idx">${i + 1}</td>
          <td class="label">${esc(r.label)}</td>
          <td class="seq mono">${esc(r.seq)}</td>
          <td colspan="7" class="err">skipped — ${esc(r.error)}</td>
        </tr>`;
      }
      const open = state.selected === i;
      return `<tr class="row${open ? " row-open" : ""}">
        <td class="idx">${i + 1}</td>
        <td class="label">${esc(r.label)}</td>
        <td class="seq mono">${esc(r.seq)}</td>
        <td>${chip(r.bbb!, "BBB")}</td>
        <td class="num">${bar(r.bbb!.prob)} <span class="p">${r.bbb!.prob.toFixed(2)}</span></td>
        <td>${chip(r.cpp!, "cell")}</td>
        <td class="num">${r.immune!.score}</td>
        <td class="num">${r.mfg!.score}</td>
        <td class="num">${pct(r.bbb!.maxSim)}</td>
        <td><button class="mini" data-open="${i}" aria-expanded="${open}">${open ? "−" : "+"}</button></td>
      </tr>
      ${open ? detailRow(r) : ""}`;
    })
    .join("");
  return `<div class="table-wrap"><table class="results">
    <thead><tr>
      <th>#</th><th>name</th><th>sequence</th><th>BBB</th><th>P(BBB)</th>
      <th>cell</th><th>immune<br><small>compat.</small></th><th>mfg</th>
      <th title="max sequence similarity to any training peptide">sim.</th><th></th>
    </tr></thead>
    <tbody>${rows}</tbody>
  </table></div>`;
}
