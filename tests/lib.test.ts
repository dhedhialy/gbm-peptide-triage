// Unit tests for the parts the UI leans on: parsing, conformal verdicts,
// similarity, dominance and the two heuristic indexes (including the
// consistency of the displayed component breakdown).

import { describe, expect, it } from "vitest";
import { parseInput, validate } from "../src/lib/parse";
import { conformalSet, verdictOf } from "../src/lib/model";
import { maxSimTo, nearestNeighbors, sim } from "../src/lib/identity";
import { paretoFront } from "../src/lib/pareto";
import { immuneCompatibility, manufacturability } from "../src/lib/heuristics";
import { ALPHAS, analyzeBatch, EXAMPLE_BATCH } from "../src/lib/triage";
import raw from "../public/model.json";
import type { ModelFile } from "../src/lib/model";

describe("input parsing", () => {
  it("reads plain sequences one per line", () => {
    const out = parseInput("TFFYGGSRGKRNNFKLKFD\nGRKKRRQRRRPPQ\n");
    expect(out.map((p) => p.seq)).toEqual(["TFFYGGSRGKRNNFKLKFD", "GRKKRRQRRRPPQ"]);
    expect(out.every((p) => !p.error)).toBe(true);
  });

  it("reads FASTA with wrapped sequence lines", () => {
    const out = parseInput(">angiopep\nTFFYGGSRGK\nRNNFKLKFD\n>x\nACDEFGHIKL");
    expect(out).toHaveLength(2);
    expect(out[0].label).toBe("angiopep");
    expect(out[0].seq).toBe("TFFYGGSRGKRNNFKLKFD");
  });

  it("lowercases are normalised, junk is rejected with a reason", () => {
    const out = parseInput("acdef\nPEP1");
    expect(out[0].seq).toBe("ACDEF");
    expect(out[0].error).toBeUndefined();
    expect(out[1].error).toMatch(/non-standard/);
  });

  it("rejects over-long input instead of scoring it", () => {
    expect(validate("A".repeat(81))).toMatch(/longer than 80/);
  });

  it("empty input yields an empty batch", () => {
    expect(parseInput("  \n\n")).toEqual([]);
  });
});

describe("conformal verdicts", () => {
  it("gives a positive call when only the positive label survives", () => {
    expect(verdictOf(conformalSet(0.95, 0.64))).toBe("positive");
  });

  it("gives a negative call when only the negative label survives", () => {
    expect(verdictOf(conformalSet(0.3, 0.64))).toBe("negative");
  });

  it("abstains when both labels survive", () => {
    expect(verdictOf(conformalSet(0.5, 0.64))).toBe("abstain");
  });

  it("abstains rather than guessing on an empty set", () => {
    expect(verdictOf(conformalSet(0.999, 0.0))).toBe("abstain");
  });

  it("a tighter quantile can only shrink the set", () => {
    for (const p of [0.05, 0.5, 0.95]) {
      expect(
        conformalSet(p, 0.1).filter(Boolean).length,
      ).toBeLessThanOrEqual(conformalSet(p, 0.6).filter(Boolean).length);
    }
  });
});

describe("sequence similarity (2·LCS/(|a|+|b|))", () => {
  const cases: [string, string, number][] = [
    ["AAAA", "AA", 0.6666666667],
    ["TFFYGGSRGKRNNFKLKFD", "TFFYGGSRGKRNNFKLKFD", 1.0],
    ["TFFYGGSRGKRNNFKLKFD", "DQKFKLNNRKGGSGYFFT", 0.3783783784],
    ["GRKKRRQRRRPPQ", "RRRRRRRRRRR", 0.5],
    ["CDIFTNSRGKRA", "CKGKGACRKTMYDCCSGSCGRRGKC", 0.3243243243],
    ["GHPGHPGARQLAGP", "GPQGPAGRLQHAG", 0.6666666667],
    ["WWWWW", "AAAAA", 0.0],
  ];
  for (const [a, b, expected] of cases) {
    it(`${a} vs ${b} = ${expected}`, () => {
      expect(sim(a, b)).toBeCloseTo(expected, 8);
      expect(sim(b, a)).toBeCloseTo(expected, 8);
    });
  }

  it("orders neighbours by similarity and reports their label", () => {
    const pool = [
      { s: "TFFYGGSRGKRNNFKLKFD", y: 1 },
      { s: "GRKKRRQRRRPPQ", y: 0 },
      { s: "WWWWWWWWWWW", y: 0 },
    ];
    const nn = nearestNeighbors("TFFYGGSRGKRNNFKLKFD", pool);
    expect(nn[0].seq).toBe("TFFYGGSRGKRNNFKLKFD");
    expect(nn[0].label).toBe(1);
    expect(nn[0].sim).toBe(1);
  });

  it("returns nothing when no training peptide shares a single residue", () => {
    expect(maxSimTo("WWWWW", [{ s: "AAAAA", y: 1 }])).toBeUndefined();
  });
});

describe("pareto front", () => {
  it("keeps only non-dominated points", () => {
    const front = paretoFront([
      { dims: [1, 1] },
      { dims: [0, 2] },
      { dims: [2, 0] },
      { dims: [0, 0] },
    ]);
    expect(front).toEqual([true, true, true, false]);
  });

  it("treats identical points as mutually non-dominated", () => {
    expect(paretoFront([{ dims: [1, 1] }, { dims: [1, 1] }])).toEqual([true, true]);
  });

  it("handles a single point and an empty batch", () => {
    expect(paretoFront([{ dims: [3, 3] }])).toEqual([true]);
    expect(paretoFront([])).toEqual([]);
  });
});

describe("heuristic indexes", () => {
  const reconcileImmune = (s: string) => {
    const h = immuneCompatibility(s);
    const sum = h.components.reduce((a, c) => a + c.points, 0);
    return { shown: 100 + sum, score: h.score };
  };
  const reconcileMfg = (s: string) => {
    const h = manufacturability(s);
    const sum = h.components.reduce((a, c) => a + c.points, 0);
    return { shown: sum, score: h.score };
  };

  const peptides = [
    "TFFYGGSRGKRNNFKLKFD",
    "GRKKRRQRRRPPQ",
    "RQIKIWFQNRRMKWKK",
    "DDDDD",
    "AAAAAAAAAAAAAAA",
    "CKGKGACRKTMYDCCSGSCGRRGKC",
    "QNGSVMWCGDGT",
  ];

  for (const s of peptides) {
    it(`immune breakdown reconciles with its score: ${s}`, () => {
      const { shown, score } = reconcileImmune(s);
      expect(score).toBe(Math.round(Math.max(0, Math.min(100, shown))));
    });

    it(`manufacturability breakdown reconciles with its score: ${s}`, () => {
      const { shown, score } = reconcileMfg(s);
      expect(score).toBe(Math.round(Math.max(0, Math.min(100, shown))));
    });
  }

  it("flags the liability motifs it claims to", () => {
    const flagged = (s: string, re: RegExp) =>
      manufacturability(s).components.some((c) => re.test(c.label) && c.points < 0);
    expect(flagged("GGNGSGGGGGGGG", /deamidation/i)).toBe(true);
    expect(flagged("GGGGGGGCGGGGG", /cysteine/i)).toBe(true);
    expect(flagged("GGGGGGGDGTGGGG", /isomerization/i)).toBe(true);
    expect(flagged("GGGGGGGGGGGGG", /Met|Trp/)).toBe(false);
    expect(flagged("QGGGGGGGGGGGGG", /N-terminal Gln/i)).toBe(true);
    const clean = manufacturability("GGGGGGGGGGGGG").components;
    expect(clean.some((c) => c.points < 0)).toBe(false);
  });

  it("scores a fully hydrophobic stretch as hard to manufacture", () => {
    expect(manufacturability("AAAAAAAAAAAAAAA").score).toBe(0);
  });

  it("scores a strongly hydrophobic 15-mer as maximally MHC-II anchored", () => {
    expect(immuneCompatibility("AAAAAAAAAAAAAAA").score).toBe(0);
  });

  it("keeps every component's points finite and labelled", () => {
    for (const s of peptides) {
      for (const h of [immuneCompatibility(s), manufacturability(s)]) {
        for (const c of h.components) {
          expect(Number.isFinite(c.points)).toBe(true);
          expect(c.label.length).toBeGreaterThan(0);
          expect(c.detail.length).toBeGreaterThan(0);
        }
      }
    }
  });
});

describe("batch triage end to end", () => {
  const model = raw as unknown as ModelFile;

  it("scores the shipped example batch at every alpha", () => {
    const peptides = parseInput(EXAMPLE_BATCH);
    expect(peptides.length).toBe(9);
    for (let ai = 0; ai < ALPHAS.length; ai++) {
      const results = analyzeBatch(model, peptides, ai);
      expect(results).toHaveLength(peptides.length);
      for (const r of results) {
        expect(r.error).toBeUndefined();
        expect(r.bbb.prob).toBeGreaterThanOrEqual(0);
        expect(r.bbb.prob).toBeLessThanOrEqual(1);
        expect(r.bbb.reasons.length).toBeGreaterThan(0);
        expect(r.bbb.neighbors.length).toBeGreaterThan(0);
        expect(r.objectives.bbb).toBeCloseTo(r.bbb.prob * 100, 9);
        expect(r.objectives.cpp).toBeCloseTo(r.cpp.prob * 100, 9);
        expect(r.objectives.immune).toBe(r.immune.score);
        expect(r.objectives.mfg).toBe(r.mfg.score);
      }
    }
  });

  it("raises coverage by abstaining more, never the other way round", () => {
    const peptides = parseInput(EXAMPLE_BATCH);
    const loose = analyzeBatch(model, peptides, 0); // alpha 0.40 -> 60% coverage
    const tight = analyzeBatch(model, peptides, ALPHAS.length - 1); // alpha 0.05 -> 95%
    const abstained = (rs: typeof loose) => rs.filter((r) => r.bbb.verdict === "abstain").length;
    expect(abstained(loose)).toBeLessThanOrEqual(abstained(tight));
  });

  it("flags length extrapolation as out-of-domain and withholds a conformal claim", () => {
    const long = "GAVLIPKGAAQTGAVLIPKGAAQTGAVLIPKGA";
    const [r] = analyzeBatch(model, parseInput(`>long\n${long}`), 4);
    expect(r.seq).toBe(long);
    expect(r.bbb.ood).toBe(true);
    expect(r.bbb.reasons.join(" ")).toMatch(/outside the 5–30 aa range/);
    expect(r.bbb.reasons.join(" ")).toMatch(/Verdict withheld/);
    expect(r.bbb.reasons.join(" ")).not.toMatch(/Conformal singleton/);
  });

  it("keeps an in-distribution peptide of the same length from being flagged", () => {
    const [r] = analyzeBatch(model, parseInput("TFFYGGSRGKRNNFKLKFD"), 4);
    expect(r.bbb.ood).toBe(false);
    expect(r.bbb.reasons.join(" ")).toMatch(/Conformal singleton|Both labels survive/);
  });

  it("reports an error entry instead of scoring an invalid peptide", () => {
    const results = analyzeBatch(model, parseInput("PEPTIDE1\nGRKKRRQRRRPPQ"), 4);
    expect(results[0].error).toMatch(/non-standard/);
    expect(results[0].bbb).toBeUndefined();
    expect(results[1].bbb).toBeDefined();
  });
});
