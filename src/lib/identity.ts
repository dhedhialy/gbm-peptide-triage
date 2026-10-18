// Sequence similarity — identical definition to scripts/train.py:
//   sim(a, b) = 2 * LCS(a, b) / (len_a + len_b)
// `floor` is the caller's decision threshold; the cheap length-bound reject
// only fires when the ratio provably cannot reach it.

function lcsLen(a: string, b: string): number {
  if (a.length < b.length) [a, b] = [b, a];
  let prev = new Array<number>(b.length + 1).fill(0);
  for (let i = 0; i < a.length; i++) {
    const cur = [0];
    for (let j = 0; j < b.length; j++) {
      cur.push(a[i] === b[j] ? prev[j] + 1 : Math.max(prev[j + 1], cur[j]));
    }
    prev = cur;
  }
  return prev[b.length];
}

export function sim(a: string, b: string, floor = 0.0): number {
  if (a === b) return 1.0;
  const lo = Math.min(a.length, b.length);
  const hi = Math.max(a.length, b.length);
  if ((2.0 * lo) / (lo + hi) < floor) return 0.0;
  return (2.0 * lcsLen(a, b)) / (a.length + b.length);
}

export interface Neighbor {
  index: number;
  seq: string;
  label: number;
  sim: number;
}

export function nearestNeighbors(
  query: string,
  pool: { s: string; y: number }[],
  k = 3,
): Neighbor[] {
  const out: Neighbor[] = [];
  for (let i = 0; i < pool.length; i++) {
    const s = sim(query, pool[i].s, 0.0);
    if (s > 0) out.push({ index: i, seq: pool[i].s, label: pool[i].y, sim: s });
  }
  out.sort((x, y) => y.sim - x.sim || x.index - y.index);
  return out.slice(0, k);
}

export function maxSimTo(query: string, pool: { s: string; y: number }[]): Neighbor | undefined {
  const nn = nearestNeighbors(query, pool, 1);
  return nn[0];
}
