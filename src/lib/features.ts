// Feature descriptors — must stay bit-for-bit identical to features() in
// scripts/train.py (verified against exported fixtures in tests).

export const AA = "ACDEFGHIKLMNPQRSTVWY";

const AA_KD: Record<string, number> = {
  A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5,
  E: -3.5, G: -0.4, H: -3.2, I: 4.5, L: 3.8, K: -3.9,
  M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9,
  Y: -1.3, V: 4.2,
};

const AA_MW: Record<string, number> = {
  A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388,
  Q: 128.1307, E: 129.1155, G: 57.0519, H: 137.1411, I: 113.1594,
  L: 113.1594, K: 128.1741, M: 131.1926, F: 147.1766, P: 97.1167,
  S: 87.0782, T: 101.1051, W: 186.2132, Y: 163.1760, V: 99.1326,
};

const PKA_ACID: Record<string, number> = { D: 3.65, E: 4.25, C: 8.33, Y: 10.1 };
const PKA_BASE: Record<string, number> = { K: 10.5, R: 12.0, H: 6.0 };
const PKA_NTERM = 9.6;
const PKA_CTERM = 2.35;
const PH = 7.4;

const HYDRO_POS = new Set(AA.split("").filter((a) => AA_KD[a] > 0));

export const BASE_FEATURES = [
  "length", "net_charge", "gravy", "aromatic_frac", "polar_frac",
  "charged_frac", "aliphatic_frac", "gly_frac", "pro_frac", "cys_frac",
  "acidic_frac", "basic_frac", "max_hydro_run", "mw",
];

export const FEATURES = BASE_FEATURES.concat(AA.split("").map((a) => `frac_${a}`));

const count = (seq: string, s: string) => {
  let n = 0;
  for (const c of seq) if (s.includes(c)) n++;
  return n;
};

export function features(seq: string): number[] {
  const n = seq.length;
  let charge = 1.0 / (1.0 + 10 ** (PH - PKA_NTERM));
  charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - PH));
  for (const aa of Object.keys(PKA_ACID)) {
    charge -= count(seq, aa) / (1.0 + 10 ** (PKA_ACID[aa] - PH));
  }
  for (const aa of Object.keys(PKA_BASE)) {
    charge += count(seq, aa) / (1.0 + 10 ** (PH - PKA_BASE[aa]));
  }
  let run = 0;
  let best = 0;
  for (const c of seq) {
    run = HYDRO_POS.has(c) ? run + 1 : 0;
    best = Math.max(best, run);
  }
  let kd = 0;
  let mw = 0;
  for (const c of seq) {
    kd += AA_KD[c];
    mw += AA_MW[c];
  }
  const f = [
    n,
    charge,
    kd / n,
    count(seq, "FWY") / n,
    count(seq, "STNQ") / n,
    count(seq, "DEKRH") / n,
    count(seq, "VLI") / n,
    count(seq, "G") / n,
    count(seq, "P") / n,
    count(seq, "C") / n,
    count(seq, "DE") / n,
    count(seq, "KRH") / n,
    best,
    mw + 18.0153,
  ];
  for (const a of AA) f.push(count(seq, a) / n);
  return f;
}

export function netCharge(seq: string): number {
  return features(seq)[1];
}
