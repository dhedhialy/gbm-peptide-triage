import { FEATURES } from "../lib/features";
import type { ModelFile } from "../lib/model";
import { pct } from "./dom";

export function renderMethods(host: HTMLElement, model: ModelFile): void {
  const gate = pct(model.ood_sim);
  host.innerHTML = `
    <section class="card">
      <h2>What this is</h2>
      <p class="lede">A batch triage tool for glioblastoma peptide candidates. It scores four
      objectives, draws the Pareto trade-offs, and — the part served by no existing BBB
      server — tells you when it does not know, with a calibrated reason.</p>
      <div class="grid-2">
        <div class="card inner">
          <h3>Trained and calibrated</h3>
          <ul class="reasons">
            <li><strong>BBB penetration</strong> — random forest on B3Pdb/B3Pred data,
            conformal prediction sets with a marginal coverage guarantee.</li>
            <li><strong>Cell penetration</strong> — same machinery on CPP-vs-random data;
            a membrane-permeability proxy, not a tumour-selectivity claim.</li>
          </ul>
        </div>
        <div class="card inner">
          <h3>Transparent heuristic indexes</h3>
          <ul class="reasons">
            <li><strong>Immune compatibility</strong> — MHC-II anchor density and class-II
            length window, with every component shown.</li>
            <li><strong>Manufacturability</strong> — charge, aggregation runs and known
            liability motifs (deamidation, oxidation, isomerisation, free cysteines).</li>
          </ul>
        </div>
      </div>
      <p class="hint">Heuristic indexes are deterministic rules, deliberately not sold as
      predictions. For calibrated MHC-II binding, run the sequences through NetMHCIIpan and
      override.</p>
    </section>

    <section class="card">
      <h2>Conformal prediction (the abstention)</h2>
      <p class="lede">Split conformal with the LAC nonconformity score. For each candidate
      sequence x and label y, score s(x, y) = 1 − p̂(y | x). On a held-out calibration set of
      n peptides, take q as the ⌈(n+1)(1−α)⌉/n empirical quantile of the true-label scores.
      The prediction set is {y : s(x, y) ≤ q}.</p>
      <ul class="reasons">
        <li>One label survives → a call is made (BBB+ / BBB−).</li>
        <li>Both survive → <strong>abstain</strong>: under the coverage target, this peptide
        cannot be separated.</li>
        <li>Coverage is marginal and in-distribution: across many draws, the true label is in
        the set at least 1−α of the time. It is a guarantee on the set, not a promise about
        any single peptide.</li>
      </ul>
      <p class="hint">The α dial in the Predict tab changes q only — the numbers next to it are
      measured on held-out clusters, not asserted. One honest caveat: this is split conformal
      on a single 158-peptide calibration set, so the coverage guarantee is marginal and
      seed-dependent, and the six-point α grid lands on discrete empirical quantiles — that is
      why the measured dial curve is not monotone step to step.</p>
    </section>

    <section class="card">
      <h2>Why the tool can refuse: out-of-domain</h2>
      <p class="lede">A conformal guarantee is void for peptides that do not resemble the
      training set. For every input we compute sequence similarity to all training peptides,
      sim(a,b) = 2·LCS(a,b) / (|a|+|b|), and compare the maximum against a gate of
      <strong>${gate}</strong>. Below that — or when the length falls outside the training
      range of 5–30 aa — the peptide is flagged out-of-domain and the conformal verdict is
      withheld: the interface says so instead of printing a confident number.</p>
      <p class="hint">The gate is not a guessed constant. It is the 10th percentile of the
      similarity from each held-out test peptide to the training pool (0.43 here, floored at
      0.40), so by construction it flags the bottom decile of peptides that already resemble
      the training data. A fixed 0.40 was tried first and measured: it flagged only 0.6% of
      held-out peptides, i.e. it never fired on realistic input. The number and its failure
      rate are both exported in <code>benchmark.json</code>.</p>
    </section>

    <section class="card">
      <h2>The models, exactly</h2>
      <p>Random forest, 200 trees, max depth 6, min 3 samples per leaf, class-balanced, on
      ${FEATURES.length} transparent descriptors per peptide:</p>
      <p class="mono wrap-text">${FEATURES.join(" · ")}</p>
      <p class="hint">No 9,000-feature Pfeature kitchen sink, no learned embeddings. The point
      of this tool is not raw accuracy — it is knowing when accuracy is not there. Full
      disclosure: our BBB forest on the similarity split reaches AUROC ≈ 0.82, below B3Pred's
      published 0.88–0.90 on random splits, and part of that gap is the split protocol itself
      (Validation tab).</p>
    </section>

    <section class="card">
      <h2>Heuristic definitions (so you can disagree precisely)</h2>
      <div class="grid-2">
        <div class="card inner">
          <h3>Immune compatibility (higher = safer)</h3>
          <p>Start at 100, subtract:</p>
          <ul class="reasons">
            <li><strong>75 ×</strong> mean MHC-II anchor density — over every 9-mer frame,
            the fraction of P1/P4/P6/P9 positions occupied by a hydrophobic residue
            (A/V/I/L/M/F/W/Y).</li>
            <li><strong>25 ×</strong> length factor — 1 inside 13–25 aa (optimal class-II
            window), 0.5 for 9–12 or 26–30, else 0.15.</li>
          </ul>
        </div>
        <div class="card inner">
          <h3>Manufacturability (higher = easier)</h3>
          <p>50 × charge score (favorable |net charge| 2–6 at pH 7.4, falling off outside)
          + 50 × aggregation score (longest hydrophobic run ≤ 4 good, ≥ 8 bad), minus
          penalties: −8 per deamidation motif (N[GSTHA]), −4 per Met, −6 per Trp, −8 per
          free cysteine, −6 per Asp isomerisation motif (D[GST]), −6 for N-terminal Gln.
          Each penalty caps at 4 occurrences.</p>
        </div>
      </div>
    </section>

    <section class="card">
      <h2>Benchmark protocol</h2>
      <ul class="reasons">
        <li><strong>Random split</strong> — stratified 60/20/20, the protocol the field
        reports.</li>
        <li><strong>Similarity split</strong> — every peptide clustered at ≥80% mutual
        LCS-similarity (union-find over all pairs); the classifier never sees a cluster in
        training that appears in test. Calibration uses separate held-out clusters.</li>
        <li><strong>B3Pred audit</strong> — our model trained on B3Pred's published train
        partitions, evaluated on their published validation partitions; then the similarity
        of every validation peptide to the training pool is measured.</li>
        <li>All numbers: 5 seeds, mean ± sd, reproduced by <code>python3 scripts/train.py</code>.</li>
      </ul>
    </section>

    <section class="card">
      <h2>Known failure cases (read before trusting anything)</h2>
      <ul class="reasons">
        <li><strong>Non-canonical residues</strong> (D-amino acids, N-methylation,
        β-amino acids, cyclisation, PEGylation) are invisible to every descriptor here — the
        input validation rejects them rather than pretending.</li>
        <li><strong>Lengths outside 5–30 aa</strong> extrapolate; the UI flags them as
        out-of-domain.</li>
        <li><strong>Scrambled sequences beat composition-only features.</strong> The example
        batch includes a scrambled Angiopep-2 control, and the model gives it a similar
        probability to the parent (0.58 at α=0.15) because 13 of the 34 descriptors are
        amino-acid composition, which a shuffle cannot change. Position-aware features are the
        fix; this tool does not have them, so treat high probabilities on shuffled or
        low-complexity sequences as uninformative.</li>
        <li><strong>Membrane-active CPPs</strong> that are not BBB shuttles can score high on
        cell penetration and moderate on BBB — judged against CPP/random negatives, not
        against the BBB endothelium in vivo.</li>
        <li><strong>Aggregation-prone leads</strong> may pass the manufacturability index:
        the hydrophobic-run rule is a screening heuristic, not a Zyggregator-class predictor.</li>
        <li><strong>Immunogenicity index ≠ immunogenicity.</strong> HLA-DRB1 alleles differ;
        a two-anchor rule cannot capture allele-specific binding.</li>
        <li><strong>Small positives (265 after removing sequences duplicated in the CPP
        negatives; 269 in the source file).</strong> The similarity-split BBB number is the one
        that matters, and it is 0.82 AUROC, not 0.95.</li>
        <li><strong>Short peptides get similarity for free.</strong> 2·LCS/(|a|+|b|) is large
        for short sequences (a 3-mer shares its whole self with any 5-mer), so the
        out-of-domain gate is weakest exactly where the feature set is also weakest.</li>
      </ul>
    </section>

    <section class="card">
      <h2>Data & citations</h2>
      <ul class="reasons">
        <li>Positive/negative data: B3Pred web server datasets (B3Pdb-derived), downloaded
        ${new Date().toISOString().slice(0, 10)}. Kumar V., Patiyal S., Dhall A., Sharma N.,
        Raghava G.P.S. <em>B3Pred: A Random-Forest-Based Method for Predicting and Designing
        Blood–Brain Barrier Penetrating Peptides.</em> Pharmaceutics 2021, 13(8):1237.</li>
        <li>Conformal prediction: Lei J. et al. <em>Distribution-Free Predictive Inference
        for Regression.</em> JASA 2018; Angelopoulos & Bates, <em>A Gentle Introduction to
        Conformal Prediction and Distribution-Free Uncertainty Quantification</em></li>
        <li>Kyte–Doolittle hydropathy (1982); LCS-based sequence identity as defined above.
        Ramírez et al. for the minimal-prefix notation this deliberately avoids.</li>
      </ul>
      <p class="hint">Not built yet, not claimed: calibrated immunogenicity (Needs
      NetMHCIIpan), design mode with single-residue edits, full-data refit. Each is a clear
      next step on top of this scaffold, none is faked here.</p>
    </section>`;
}
