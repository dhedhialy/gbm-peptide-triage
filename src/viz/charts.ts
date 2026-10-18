// Minimal hand-rolled SVG charts. Monochrome by construction: shape, fill and
// dash encode state; interaction (hover/click/crosshair) carries the detail.

const NS = "http://www.w3.org/2000/svg";

function el<K extends keyof SVGElementTagNameMap>(
  tag: K,
  attrs: Record<string, string | number> = {},
  parent?: Element,
): SVGElementTagNameMap[K] {
  const node = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, String(v));
  parent?.appendChild(node);
  return node;
}

export interface Tooltip {
  show(x: number, y: number, html: string): void;
  hide(): void;
}

export function makeTooltip(host: HTMLElement): Tooltip {
  const div = document.createElement("div");
  div.className = "chart-tooltip";
  host.appendChild(div);
  return {
    show(x, y, html) {
      div.innerHTML = html;
      div.style.opacity = "1";
      div.style.left = `${x}px`;
      div.style.top = `${y}px`;
    },
    hide() {
      div.style.opacity = "0";
    },
  };
}

interface Margin {
  l: number;
  r: number;
  t: number;
  b: number;
}

function frame(host: HTMLElement, width: number, height: number, m: Margin) {
  const svg = el("svg", {
    viewBox: `0 0 ${width} ${height}`,
    class: "chart",
    style: `max-width:${width}px`,
  });
  host.appendChild(svg);
  return {
    svg,
    x: (v: number, dom: [number, number]) =>
      m.l + ((v - dom[0]) / (dom[1] - dom[0] || 1)) * (width - m.l - m.r),
    y: (v: number, dom: [number, number]) =>
      height - m.b - ((v - dom[0]) / (dom[1] - dom[0] || 1)) * (height - m.t - m.b),
  };
}

function axes(
  svg: SVGElementTagNameMap["svg"],
  m: Margin,
  width: number,
  height: number,
  xDom: [number, number],
  yDom: [number, number],
  xLabel: string,
  yLabel: string,
  xTickFmt: (v: number) => string = (v) => String(Math.round(v)),
) {
  const px = (v: number) => m.l + ((v - xDom[0]) / (xDom[1] - xDom[0] || 1)) * (width - m.l - m.r);
  const py = (v: number) => height - m.b - ((v - yDom[0]) / (yDom[1] - yDom[0] || 1)) * (height - m.t - m.b);
  for (let i = 0; i <= 4; i++) {
    const xv = xDom[0] + ((xDom[1] - xDom[0]) * i) / 4;
    const yv = yDom[0] + ((yDom[1] - yDom[0]) * i) / 4;
    el("line", { x1: m.l, x2: width - m.r, y1: py(yv), y2: py(yv), class: "grid" }, svg);
    el("line", { x1: px(xv), x2: px(xv), y1: m.t, y2: height - m.b, class: "grid" }, svg);
    el("text", { x: px(xv), y: height - m.b + 16, class: "tick", "text-anchor": "middle" }, svg).textContent = xTickFmt(xv);
    el("text", { x: m.l - 8, y: py(yv) + 3, class: "tick", "text-anchor": "end" }, svg).textContent = xTickFmt(yv);
  }
  el("line", { x1: m.l, x2: width - m.r, y1: height - m.b, y2: height - m.b, class: "axis" }, svg);
  el("line", { x1: m.l, x2: m.l, y1: m.t, y2: height - m.b, class: "axis" }, svg);
  el("text", { x: (m.l + width - m.r) / 2, y: height - 6, class: "axis-label", "text-anchor": "middle" }, svg).textContent = xLabel;
  const yl = el(
    "text",
    { x: 12, y: (m.t + height - m.b) / 2, class: "axis-label", "text-anchor": "middle", transform: `rotate(-90 12 ${(m.t + height - m.b) / 2})` },
    svg,
  );
  yl.textContent = yLabel;
}

// ------------------------------------------------------------------ scatter --

export interface ScatterPoint {
  x: number;
  y: number;
  label: string;
  sub?: string;
  kind: "front" | "plain" | "abstain" | "ghost";
}

export interface ScatterCfg {
  points: ScatterPoint[];
  xLabel: string;
  yLabel: string;
  xDomain: [number, number];
  yDomain: [number, number];
  picked?: number;
  onPick?: (i: number) => void;
  width?: number;
  height?: number;
}

export function scatter(host: HTMLElement, cfg: ScatterCfg): void {
  host.querySelector("svg")?.remove();
  const width = cfg.width ?? 640;
  const height = cfg.height ?? 420;
  const m: Margin = { l: 46, r: 16, t: 14, b: 40 };
  const { svg, x, y } = frame(host, width, height, m);
  axes(svg, m, width, height, cfg.xDomain, cfg.yDomain, cfg.xLabel, cfg.yLabel);
  const tip = makeTooltip(host);

  const groups: SVGElementTagNameMap["g"][] = [];
  cfg.points.forEach((p, i) => {
    if (p.kind === "ghost") {
      el("circle", { cx: x(p.x, cfg.xDomain), cy: y(p.y, cfg.yDomain), r: 2.2, class: "pt-ghost" }, svg);
    } else {
      const g = el("g", { class: `pt pt-${p.kind}`, tabindex: 0 }, svg);
      const cx = x(p.x, cfg.xDomain);
      const cy = y(p.y, cfg.yDomain);
      el("circle", { cx, cy, r: p.kind === "front" ? 5.5 : 4.5 }, g);
      if (i === cfg.picked) el("circle", { cx, cy, r: 9.5, class: "pt-ring" }, g);
      groups.push(g);
      const enter = () => {
        tip.show(
          cx + 12,
          cy - 8,
          `<strong>${p.label}</strong>${p.sub ? `<br>${p.sub}` : ""}`,
        );
      };
      g.addEventListener("mouseenter", enter);
      g.addEventListener("focus", enter);
      g.addEventListener("mouseleave", () => tip.hide());
      if (cfg.onPick) {
        g.addEventListener("click", () => cfg.onPick!(i));
        g.addEventListener("keydown", (e) => {
          if (e.key === "Enter" || e.key === " ") cfg.onPick!(i);
        });
      }
    }
  });
  // nearest-point follow: hovering anywhere shows the closest peptide
  const rect = el("rect", { x: m.l, y: m.t, width: width - m.l - m.r, height: height - m.t - m.b, fill: "transparent" }, svg);
  rect.addEventListener("mousemove", (ev) => {
    const svgRect = svg.getBoundingClientRect();
    const sx = ((ev as MouseEvent).clientX - svgRect.left) * (width / svgRect.width);
    const sy = ((ev as MouseEvent).clientY - svgRect.top) * (height / svgRect.height);
    let best = -1;
    let bestD = 30;
    cfg.points.forEach((p, i) => {
      if (p.kind === "ghost") return;
      const d = Math.hypot(x(p.x, cfg.xDomain) - sx, y(p.y, cfg.yDomain) - sy);
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    });
    if (best >= 0) {
      const p = cfg.points[best];
      tip.show(
        x(p.x, cfg.xDomain) + 12,
        y(p.y, cfg.yDomain) - 8,
        `<strong>${p.label}</strong>${p.sub ? `<br>${p.sub}` : ""}`,
      );
    } else tip.hide();
  });
  rect.addEventListener("mouseleave", () => tip.hide());
}

// --------------------------------------------------------------- line chart --

export interface LineSeries {
  name: string;
  points: { x: number; y: number }[];
  dash?: boolean;
  markers?: boolean;
}

export interface LineCfg {
  series: LineSeries[];
  xLabel: string;
  yLabel: string;
  xDomain: [number, number];
  yDomain: [number, number];
  diagonal?: boolean;
  width?: number;
  height?: number;
}

export function lineChart(host: HTMLElement, cfg: LineCfg): void {
  host.querySelector("svg")?.remove();
  const width = cfg.width ?? 420;
  const height = cfg.height ?? 320;
  const m: Margin = { l: 46, r: 16, t: 14, b: 40 };
  const { svg, x, y } = frame(host, width, height, m);
  axes(svg, m, width, height, cfg.xDomain, cfg.yDomain, cfg.xLabel, cfg.yLabel, (v) => v.toFixed(2));
  const tip = makeTooltip(host);
  if (cfg.diagonal) {
    el("line", {
      x1: x(cfg.xDomain[0], cfg.xDomain), y1: y(cfg.yDomain[0], cfg.yDomain),
      x2: x(cfg.xDomain[1], cfg.xDomain), y2: y(cfg.yDomain[1], cfg.yDomain),
      class: "diag",
    }, svg);
  }
  cfg.series.forEach((s, si) => {
    const d = s.points
      .map((p, i) => `${i === 0 ? "M" : "L"} ${x(p.x, cfg.xDomain).toFixed(1)} ${y(p.y, cfg.yDomain).toFixed(1)}`)
      .join(" ");
    el("path", { d, class: `series series-${si}${s.dash ? " series-dash" : ""}` }, svg);
    if (s.markers !== false) {
      s.points.forEach((p, i) => {
        el("circle", {
          cx: x(p.x, cfg.xDomain), cy: y(p.y, cfg.yDomain),
          r: i === s.points.length - 1 ? 3.8 : 3.2, class: "marker",
        }, svg);
      });
    }
  });
  const cross = el("line", { class: "crosshair", y1: m.t, y2: height - m.b, x1: -10, x2: -10 }, svg);
  const rect = el("rect", { x: m.l, y: m.t, width: width - m.l - m.r, height: height - m.t - m.b, fill: "transparent" }, svg);
  rect.addEventListener("mousemove", (ev) => {
    const svgRect = svg.getBoundingClientRect();
    const sx = ((ev as MouseEvent).clientX - svgRect.left) * (width / svgRect.width);
    const xv = cfg.xDomain[0] + ((sx - m.l) / (width - m.l - m.r)) * (cfg.xDomain[1] - cfg.xDomain[0]);
    cross.setAttribute("x1", String(sx));
    cross.setAttribute("x2", String(sx));
    const lines = cfg.series
      .map((s) => {
        let nearest = s.points[0];
        for (const p of s.points) if (Math.abs(p.x - xv) < Math.abs(nearest.x - xv)) nearest = p;
        return `${s.name}: <strong>${nearest.y.toFixed(3)}</strong> at ${nearest.x.toFixed(2)}`;
      })
      .join("<br>");
    tip.show(sx + 12, m.t + 6, lines);
  });
  rect.addEventListener("mouseleave", () => {
    cross.setAttribute("x1", "-10");
    cross.setAttribute("x2", "-10");
    tip.hide();
  });
  const legend = document.createElement("div");
  legend.className = "legend";
  legend.innerHTML = cfg.series
    .map((s) => `<span><i class="swatch${s.dash ? " swatch-dash" : ""}"></i>${s.name}</span>`)
    .join("");
  host.appendChild(legend);
}

// ------------------------------------------------------------ grouped bars --

export interface BarGroup {
  label: string;
  values: number[];
}

export interface BarsCfg {
  groups: BarGroup[];
  seriesNames: string[];
  max: number;
  yLabel: string;
  width?: number;
  height?: number;
  annotations?: string[];
}

export function groupedBars(host: HTMLElement, cfg: BarsCfg): void {
  host.querySelector("svg")?.remove();
  const width = cfg.width ?? 640;
  const height = cfg.height ?? 320;
  const m: Margin = { l: 46, r: 16, t: 14, b: 62 };
  const { svg } = frame(host, width, height, m);
  axes(svg, m, width, height, [0, cfg.groups.length], [0, cfg.max], "", cfg.yLabel, (v) =>
    Number.isInteger(v) ? String(v) : v.toFixed(1),
  );
  const tip = makeTooltip(host);
  const gw = (width - m.l - m.r) / cfg.groups.length;
  const n = cfg.seriesNames.length;
  const bw = Math.min(34, (gw - 14) / n);
  cfg.groups.forEach((g, gi) => {
    const x0 = m.l + gi * gw + (gw - bw * n) / 2;
    g.values.forEach((v, si) => {
      const h = ((height - m.t - m.b) * v) / cfg.max;
      const x = x0 + si * bw;
      const y = height - m.b - h;
      const rect = el(
        "rect",
        { x, y, width: bw - 3, height: h, class: `bar bar-${si}` },
        svg,
      );
      const enter = () => {
        tip.show(
          x + bw / 2 + 10,
          y - 6,
          `<strong>${g.label}</strong><br>${cfg.seriesNames[si]}: ${v.toFixed(3)}`,
        );
      };
      rect.addEventListener("mouseenter", enter);
      rect.addEventListener("mouseleave", () => tip.hide());
    });
    el("text", { x: m.l + gi * gw + gw / 2, y: height - m.b + 16, class: "tick", "text-anchor": "middle" }, svg).textContent = g.label;
    if (cfg.annotations?.[gi]) {
      el(
        "text",
        { x: m.l + gi * gw + gw / 2, y: height - m.b + 32, class: "tick tick-note", "text-anchor": "middle" },
        svg,
      ).textContent = cfg.annotations[gi];
    }
  });
  const legend = document.createElement("div");
  legend.className = "legend";
  legend.innerHTML = cfg.seriesNames
    .map((s, i) => `<span><i class="swatch swatch-${i}"></i>${s}</span>`)
    .join("");
  host.appendChild(legend);
}
