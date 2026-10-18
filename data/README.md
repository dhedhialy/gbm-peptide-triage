# Data provenance

All FASTA files in `data/raw/` were downloaded from the **B3Pred** web server
(<https://webs.iiitd.edu.in/raghava/b3pred/>) on 2026-10-17, which mirrors the
datasets from:

> Kumar, V., Patiyal, S., Dhall, A., Sharma, N. & Raghava, G.P.S. (2021).
> *B3Pred: A Random-Forest-Based Method for Predicting and Designing Blood–Brain
> Barrier Penetrating Peptides.* Pharmaceutics 13(8):1237.

| File | Records | Role |
|---|---|---|
| `Positive_B3PPs.fasta` | 269 | Experimentally validated BBB-penetrating peptides |
| `pos_train.fasta` / `pos_valid.fasta` | 215 / 54 | B3Pred's own split of the positives |
| `Negative_CPPs.fasta` | 269 | Cell-penetrating peptides that do **not** cross the BBB (hard negatives, set D1) |
| `Negative_Random_Balanced.fasta` | 269 | Random negatives, size-matched (set D2) |
| `Negative_Random.fasta` | 2690 | Full random negatives (set D3) |
| `*_neg_train.fasta` / `*_neg_valid.fasta` | — | B3Pred's own train/valid split for D1/D2/D3 |

License/attribution: dataset redistribution follows B3Pred's terms of use on
their download page; cite the B3Pred paper when using this data.

The BBB positives are small (269 peptides). This is the central limitation of
the whole tool and is stated on the app's Validation tab and in the README.
