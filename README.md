# Peptide triage for glioblastoma — with a conformal *I don't know*

**Live: <https://dhedhialy.github.io/gbm-peptide-triage/>**

A browser-side triage tool for blood–brain-barrier (BBB) and cell-penetrating
peptide candidates. It scores four objectives, draws the Pareto trade-offs
between them, and does the thing no BBB web server does: **it refuses to answer
when it should**, with a calibrated reason.

Everything runs client-side. No sequence you paste leaves the page.

```bash
nvm use            # Node 22+ (see package.json)
yarn install
yarn dev           # http://localhost:5173
```

The shipped `public/model.json` is already trained, so the app works on first
run. `yarn train` (Python 3, scikit-learn, numpy) reproduces it from
`data/raw/*.fasta` in ~75 s and rewrites `model.json` + `benchmark.json`.

`bash scripts/deploy_pages.sh` rebuilds and pushes `dist/` to the `gh-pages`
branch that GitHub Pages serves — branch-based deploy, so no Actions workflow
scope on the clone's token is needed.

## What it actually does

| Tab | What is there |
|---|---|
| **Predict** | Batch triage table. Each peptide gets a BBB verdict, a cell-penetration verdict, an immune-compatibility index, a manufacturability index, its nearest known peptides, and a per-feature breakdown of *why* the forest moved its probability. Rows expand. |
| **Pareto** | The batch on any two of the four objectives, with the non-dominated front, a reference cloud of training BBB+ peptides, and explicit dominator lists. Out-of-domain peptides are plotted but excluded from ranking. |
| **Validation** | Random split vs 80%-similarity cluster split; an audit of B3Pred's own published splits showing how much of their external test set is a near-duplicate of their training set; measured abstention/error curves; conformal accounting. |
| **Methods** | Exact descriptors, the conformal recipe, the out-of-domain gate and how it was derived, and a list of known failure cases. |

## The honesty core

1. **Conformal prediction, not a threshold.** Split conformal with the LAC
   score (`s = 1 − p̂(y|x)`), `q` = empirical quantile of true-label scores on a
   held-out calibration set. Prediction set `{y : s(x,y) ≤ q}`. Two survivors
   means the tool abstains. The α dial moves `q` only, and the numbers beside it
   are measured on held-out clusters.
2. **An out-of-domain gate that actually fires.** The gate is the 10th
   percentile of held-out similarity to the training pool (0.43 here, floored at
   0.40) *plus* any length outside the training range of 5–30 aa. A fixed 0.40
   was tried first and measured: it flagged 0.6% of held-out peptides, i.e. it
   never fired on realistic input. Flagged peptides get no conformal claim and
   are excluded from the Pareto ranking.
3. **The split is the result.** B3Pred's published external-validation accuracy
   is computed on splits where 47.7% / 25.0% / 5.1% of the validation peptides
   (D1/D2/D3) are ≥80% similar to a training peptide. Our AUROC on their exact
   splits is 0.900 / 0.879 / 0.916; restricted to validation peptides below 80%
   similarity it drops to 0.843 / 0.829 / 0.880. We could not re-run B3Pred
   itself under novel-only splits — its standalone package references a
   server-side model file that is not distributed — so their numbers are quoted,
   not reproduced.
4. **Heuristics are labelled heuristics.** Immune compatibility and
   manufacturability are deterministic indexes over published rules, shown with
   every component and its points. They are not calibrated models and the UI
   says so.

## Measured results

All numbers come from `scripts/train.py`, 5 seeds, mean ± sd, exported to
`public/benchmark.json` and rendered on the Validation tab.

| Task | Split | AUROC | AUPRC | Conformal coverage | Abstained | Accuracy when answering |
|---|---|---|---|---|---|---|
| BBB penetration | random 60/20/20 | 0.838 ± 0.033 | 0.754 | — | — | — |
| BBB penetration | 80%-similarity clusters | 0.815 ± 0.027 | 0.730 | 89.8% | 39.5% | 83.2% |
| Cell penetration | random 60/20/20 | 0.932 ± 0.019 | 0.944 | — | — | — |
| Cell penetration | 80%-similarity clusters | 0.924 ± 0.025 | 0.931 | 94.8% | 25.8% | 93.1% |

Dataset: 269 BBB-positive peptides, 269 CPP negatives, 269 size-matched random
negatives (265 / 256 / 269 after removing sequences that appear with conflicting
labels across files). Lengths 5–30 aa. See `data/README.md` for provenance.

## Known failure cases

The app lists these on its Methods tab; the short version: composition-only
descriptors cannot tell a scrambled peptide from its parent (the example batch
includes a scrambled Angiopep-2 control and the model scores it like the
parent), non-canonical residues are rejected rather than silently mangled,
strongly hydrophobic stretches score badly on manufacturability, the immune
index is a two-anchor rule and not HLA-allele-specific, and every metric rests
on 265 positives, so the similarity-split number (0.82 AUROC) is the one that
matters.

## Layout

```
scripts/train.py        features, splits, conformal, forest export, and the B3Pred audit
public/model.json       exported forest + quantiles + training sequences + parity fixtures
public/benchmark.json   every measured number the UI displays
src/lib/features.ts     34 descriptors — mirrors features() in train.py
src/lib/identity.ts     2·LCS/(|a|+|b|) similarity, neighbour search, OOD gate
src/lib/model.ts        forest inference, conformal sets, attribution, verdict reasons
src/lib/heuristics.ts   immune-compatibility + manufacturability indexes
src/lib/triage.ts       batch scoring and the four objectives
src/lib/pareto.ts       non-dominated front
src/viz/charts.ts       hand-rolled SVG scatter / line / grouped bars (monochrome)
src/ui/*.ts             one renderer per tab
tests/                  114 tests: Python↔TypeScript parity, library units, UI and boot
```

`yarn test` runs everything, including parity tests that re-score the exported
Python fixtures through the TypeScript implementation (features to 1e-7,
probabilities and attribution to 1e-5).

## Design note

The interface is deliberately monochrome. State is encoded by fill, outline,
dash and position rather than hue, so the figures survive greyscale printing and
read the same for colour-blind readers.

## What is not built

Calibrated immunogenicity (needs NetMHCIIpan), a design mode with
single-residue edits, and a full-data refit. Each is a clear next step on top of
this scaffold; none is faked.

## Data & citations

- Positives/negatives: B3Pred web server datasets (B3Pdb-derived).
  Kumar V., Patiyal S., Dhall A., Sharma N., Raghava G.P.S. *B3Pred: A
  Random-Forest-Based Method for Predicting and Designing Blood–Brain Barrier
  Penetrating Peptides.* Pharmaceutics 2021, 13(8):1237.
- Conformal prediction: Lei J. et al., *Distribution-Free Predictive Inference
  for Regression*, JASA 2018; Angelopoulos & Bates, *A Gentle Introduction to
  Conformal Prediction and Distribution-Free Uncertainty Quantification*.
- Kyte–Doolittle hydropathy (1982).

## License

MIT (see `LICENSE`). The datasets keep B3Pred's terms of use — cite the paper
if you use them.

Not a clinical tool. It is a screening aid that shows its reasoning and refuses
when it should.
