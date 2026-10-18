// Parity tests: the browser must reproduce the exported Python model exactly.
// Every fixture in public/model.json was produced by scripts/train.py and
// re-scored here through src/lib/*.ts.

import { describe, expect, it } from "vitest";
import raw from "../public/model.json";
import { features } from "../src/lib/features";
import { attribution, conformalSet, forestProb, type ModelFile } from "../src/lib/model";
import { maxSimTo } from "../src/lib/identity";

const model = raw as unknown as ModelFile;

describe("feature parity with train.py", () => {
  for (const fx of model.fixtures) {
    it(`${fx.task} features: ${fx.s.slice(0, 12)}…`, () => {
      const x = features(fx.s);
      expect(x.length).toBe(model.features.length);
      x.forEach((v, i) => expect(v).toBeCloseTo(fx.features[i], 7));
    });
  }
});

describe("forest inference parity", () => {
  for (const fx of model.fixtures) {
    it(`${fx.task} P(y=1) for ${fx.s.slice(0, 12)}…`, () => {
      const task = model.tasks[fx.task as "bbb" | "cpp"];
      const x = features(fx.s);
      expect(forestProb(task, x)).toBeCloseTo(fx.prob, 5);
    });
  }

  it("probabilities leave the tree value range untouched for out-of-range inputs", () => {
    const task = model.tasks.bbb;
    const p = forestProb(task, new Array(model.features.length).fill(-1e9));
    expect(p).toBeGreaterThanOrEqual(0);
    expect(p).toBeLessThanOrEqual(1);
    const q = forestProb(task, new Array(model.features.length).fill(1e9));
    expect(q).toBeGreaterThanOrEqual(0);
    expect(q).toBeLessThanOrEqual(1);
  });
});

describe("conformal set parity", () => {
  for (const fx of model.fixtures) {
    it(`${fx.task} set at the shipped alpha`, () => {
      const task = model.tasks[fx.task as "bbb" | "cpp"];
      expect(conformalSet(fx.prob, task.q)).toEqual(fx.set);
    });
  }
});

describe("similarity parity", () => {
  for (const fx of model.fixtures) {
    it(`${fx.task} max similarity to train: ${fx.s.slice(0, 12)}…`, () => {
      const task = model.tasks[fx.task as "bbb" | "cpp"];
      const nn = maxSimTo(fx.s, task.train_seqs);
      expect(nn?.sim ?? 0).toBeCloseTo(fx.sim_to_train, 4);
    });
  }
});

describe("deployed attribution matches the export", () => {
  for (const fx of model.fixtures) {
    it(`${fx.task} top attribution deltas`, () => {
      const task = model.tasks[fx.task as "bbb" | "cpp"];
      const got = attribution(task, features(fx.s)).slice(0, fx.top_attr.length);
      expect(got.length).toBe(fx.top_attr.length);
      fx.top_attr.forEach(([featureIndex, delta], k) => {
        expect(got[k][0]).toBe(featureIndex);
        expect(got[k][1]).toBeCloseTo(delta, 5);
      });
    });
  }
});

describe("alpha dial", () => {
  it("a higher coverage target needs a larger quantile, so sets only grow", () => {
    const task = model.tasks.bbb;
    const prob = 0.42;
    const loose = task.qs["0.4"]; // 60% coverage target
    const tight = task.qs["0.05"]; // 95% coverage target
    expect(tight).toBeGreaterThanOrEqual(loose);
    const size = (q: number) => conformalSet(prob, q).filter(Boolean).length;
    expect(size(tight)).toBeGreaterThanOrEqual(size(loose));
    expect(Object.keys(task.qs).length).toBe(6);
  });
});
