// Pareto dominance over maximization objectives (all four objectives are
// oriented higher-is-better before this is called).
// ponytail: O(n^2) front; batch inputs are tens of peptides, not millions.

export interface ParetoPoint {
  dims: number[];
}

export function paretoFront(points: ParetoPoint[]): boolean[] {
  const n = points.length;
  const onFront = new Array<boolean>(n).fill(true);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i === j || !onFront[i]) continue;
      if (dominates(points[j].dims, points[i].dims)) onFront[i] = false;
    }
  }
  return onFront;
}

function dominates(a: number[], b: number[]): boolean {
  let strictlyBetter = false;
  for (let k = 0; k < a.length; k++) {
    if (a[k] < b[k] - 1e-9) return false;
    if (a[k] > b[k] + 1e-9) strictlyBetter = true;
  }
  return strictlyBetter;
}
