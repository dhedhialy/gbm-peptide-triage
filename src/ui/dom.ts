export function esc(s: string): string {
  return s.replace(/[&<>"']/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]!,
  );
}

export function delegate<T extends Element>(
  root: ParentNode,
  selector: string,
  event: string,
  fn: (el: T, ev: Event) => void,
): void {
  root.querySelectorAll(selector).forEach((el) => {
    el.addEventListener(event, (ev) => fn(el as T, ev));
  });
}

export function pct(v: number, digits = 0): string {
  return `${(v * 100).toFixed(digits)}%`;
}
