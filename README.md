# Summary — `eeg_ds001787_psd_final`

**Author:** Giulia Finotto — Department of bioengineering, University of Pavia
**Dataset:** OpenNeuro ds001787 — Meditation EEG
**Main library:** MNE-Python v1.11.0 / MNE-BIDS v0.17.0

---

## Overview

This notebook implements a complete EEG preprocessing and spectral analysis pipeline for the OpenNeuro dataset **ds001787**, a multi-session meditation study comprising 24 participants (12 expert meditators, 12 novices) recorded with a 64-channel BioSemi ActiveTwo system. The pipeline follows the BIDS standard throughout and is fully reproducible: all tunable parameters are centralised in a single configuration cell.

---

## 1. Data Acquisition

The dataset is downloaded from [openneuro.org](https://openneuro.org/datasets/ds001787) using the official `openneuro-py` client. Files already present locally are skipped automatically (incremental download). The directory structure follows BIDS v1.6:

```
ds001787/
├── participants.tsv          ← one row per subject (age, gender, group)
├── dataset_description.json
└── sub-001/
    ├── ses-01/eeg/sub-001_ses-01_task-meditation_eeg.bdf
    └── ses-02/eeg/sub-001_ses-02_task-meditation_eeg.bdf
```

All 24 subject IDs are discovered automatically by reading `participants.tsv`. Each subject has up to three sessions (`ses-01`, `ses-02`, `ses-03`); missing sessions are skipped with a warning rather than aborting the pipeline.

---

## 2. Global Parameters

All tunable thresholds are defined in a single cell and propagate automatically to every downstream function. Key values:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `L_FREQ` / `H_FREQ` | 1 – 40 Hz | Band-pass filter range |
| `FLAT_STD_THRESHOLD` | 0.5 µV | Flat channel detection (std) |
| `FLAT_PTP_THRESHOLD` | 1.0 µV | Flat channel detection (ptp) |
| `EPOCH_DURATION` | 2 s | Artifact rejection window |
| `Z_THRESHOLD` | 3.0 | Artifact rejection cutoff |
| `N_FFT` | 2048 | Welch window length |
| `N_OVERLAP` | 512 | Welch window overlap (25%) |

---

## 3. Preprocessing Pipeline

Six sequential functions are defined (sections 5a–5f) and then applied to every subject in a main loop (section 6). Each function has a single responsibility, making the pipeline modular and easy to extend.

### Step 1 — Data Loading and Multi-session Concatenation (`load_subject`)

BIDS files are read with `mne_bids.read_raw_bids()`, which memory-maps the signal without loading it into RAM. Where multiple sessions are available, raw objects are concatenated along the time axis:

$$X_{\text{concat}} = \left[X^{(1)} \mid X^{(2)} \mid \cdots \mid X^{(K)}\right] \in \mathbb{R}^{C \times \sum_k T_k}$$

Session boundaries are marked as `BAD boundary` annotations so that downstream steps never straddle them.

### Step 2 — Non-EEG Channel Removal (`exclude_non_eeg_chs`)

Fifteen physiological channels (GSR1, GSR2, Resp, Erg1, Erg2, EXG1–EXG8, Plet, Temp) are dropped before re-referencing. Including them in the average reference computation would contaminate the reference mean $\bar{x}[n]$ with non-neural variance, biasing every EEG channel.

### Step 3 — Flat Channel Detection (`remove_flat_chs`)

Channels are flagged as disconnected if **both** conditions hold simultaneously across the full recording:

$$\sigma_c < 0.5\,\mu\text{V} \quad \text{AND} \quad \text{ptp}_c < 1.0\,\mu\text{V}$$

The dual criterion prevents channels with isolated transient spikes (high ptp, low std) from being misclassified as flat. Flagged channels are physically removed from the data matrix. No flat channels were detected in any of the 24 subjects.

### Step 4 — Average Re-referencing

The instantaneous scalp mean is subtracted from every channel, making the recording independent of the physical reference electrode:

$$\tilde{X} = \left(I - \frac{1}{C}\mathbf{1}\mathbf{1}^\top\right)X$$

The centering matrix $\left(I - \frac{1}{C}\mathbf{1}\mathbf{1}^\top\right)$ projects the data onto the null space of $\mathbf{1}^\top$. Re-referencing is applied twice: once before filtering (on the raw signal) and once after artifact rejection (on the clean signal), so that discarded epochs do not distort the reference mean.

### Step 5 — Band-pass FIR Filtering (`preprocess`)

A linear-phase FIR filter with a **Blackman window** retains frequencies in $[1, 40]$ Hz:

$$Y(f) = H(f) \cdot X(f), \qquad H(f) = \begin{cases} 1 & 1 \leq f \leq 40\,\text{Hz} \\ 0 & \text{otherwise} \end{cases}$$

- **High-pass at 1 Hz** removes sub-hertz baseline drift (electrode movement, sweat, DC offset).
- **Low-pass at 40 Hz** removes broadband muscular EMG artefacts and attenuates 50 Hz power-line interference (the Blackman window provides −74 dB stopband attenuation, eliminating Gibbs ripple).
- The linear-phase property guarantees constant group delay: all frequency components are delayed by the same number of samples, preserving ERP morphology.

### Step 6 — Artifact Rejection (`reject_artifacts`)

Large transient artefacts (eye blinks, jaw clenches, electrode pops) are removed using per-channel z-score thresholding on peak-to-peak amplitude:

1. The signal is divided into $N$ non-overlapping 2-second epochs of $L = f_s \times 2$ samples.
2. Peak-to-peak amplitude is computed per epoch per channel: $p_{c,k} = \max_t x_c^{(k)} - \min_t x_c^{(k)}$
3. Each channel's ptp values are z-scored independently: $z_{c,k} = (p_{c,k} - \mu_c) / \sigma_c$
4. An epoch is rejected if **any** channel exceeds $z > 3$ (logical OR across channels).
5. Surviving epochs are concatenated into a clean `RawArray`.

Per-channel normalisation removes the effect of naturally different electrode amplitudes across the scalp. The average rejection rate across subjects was **14.3%** (range: 5.4–27.6%); no subject exceeded the 40% high-concern threshold.

---

## 4. Power Spectral Density Estimation

PSD is estimated on the fully preprocessed signal using **Welch's method**:

$$\hat{S}(f_k) = \frac{1}{M}\sum_{m=0}^{M-1} P_m(f_k), \qquad \Delta f = \frac{f_s}{L} \approx 0.125\,\text{Hz}$$

Each segment is Hann-windowed to reduce spectral leakage; averaging $M$ periodograms reduces the standard error by $1/\sqrt{M}$. Mean band power for each canonical band is computed as:

$$\bar{P}_{\text{band}} = \frac{1}{C|\mathcal{F}|}\sum_{c=1}^{C}\sum_{f_k \in \mathcal{F}} \hat{S}_c(f_k)$$

and converted to decibels: $P_{\text{dB}} = 10\log_{10}(\bar{P}_{\text{band}})$.

| Band | Range (Hz) | Functional correlates |
|------|-----------|----------------------|
| Delta (δ) | 0.5 – 4 | Deep sleep, large slow waves |
| Theta (θ) | 4 – 8 | Drowsiness, early meditation, memory encoding |
| Alpha (α) | 8 – 13 | Relaxed wakefulness, meditative attentional gating |
| Beta (β) | 13 – 30 | Active cognition, focused attention |
| Gamma (γ) | 30 – 40 | High-level sensory binding, cross-region synchrony |

---

## 5. Group-level Outputs

The main loop (section 6) processes all 24 subjects successfully (0 failures). Results are accumulated in a DataFrame and exported to `results_all_subjects.csv`. The following outputs are produced:

| Output file | Content |
|-------------|---------|
| `results_all_subjects.csv` | One row per subject: channels, flat count, rejection rate, band powers |
| `psd_raw_per_subject.png` | Raw PSD panel grid — one subplot per subject (QC visualization) |
| `psd_before_after_grandavg.png` | Grand-average PSD before vs after preprocessing (±1 SD shading) |
| `psd_all_subjects.png` | Individual clean PSD curves overlaid with grand mean |
| `band_power_all_subjects.png` | Grouped bar chart of band power per subject |
| `band_power_boxplot.png` | Distribution of band power across subjects per band |
| `qc_rejection_rate.png` | Artifact rejection rate per subject (colour-coded by severity) |

---

## 6. Quality Control

Visual inspection of the raw PSD panel (section 8) identified **sub-007** as anomalous: its spectral shape was inconsistent with physiological EEG. This subject was excluded from all subsequent analyses, leaving a final analytic sample of **N = 23**.

The artifact rejection QC chart (section 13) colour-codes each subject:
- 🟢 **Green** (< 20% rejected): 17 subjects — acceptable data quality
- 🟠 **Orange** (20–40% rejected): 6 subjects — moderate; flagged for sensitivity analysis
- 🔴 **Red** (> 40% rejected): 0 subjects



# Summary — `eeg_ds001787_psd_statistical_analysis.ipynb`

**Author:** Giulia Finotto — Department of Bioengineering, University of Pavia
**Dataset:** OpenNeuro ds001787 — Meditation EEG
**Prerequisites:** `results_all_subjects.csv` (produced by `eeg_ds001787_psd_final`) · `ds001787/participants.tsv`

---

## Overview

This notebook performs the full statistical analysis of EEG band power differences between expert and novice meditators. It takes as input the per-subject band power values computed by the preprocessing pipeline and produces a publication-ready results table, three figures, and an exported CSV. All statistical procedures are implemented from scratch using only NumPy and SciPy, with no dependency on external statistics packages.

---

## Notebook Structure

| Section | Content |
|---------|---------|
| 0 | Library imports |
| 1 | Data loading: `results_all_subjects.csv` + `participants.tsv` |
| 2 | Merge, group assignment, and QC exclusion |
| 3 | Descriptive statistics by group |
| 4 | Welch's t-test for each frequency band |
| 5 | Benjamini–Hochberg FDR correction |
| 6 | Effect size: Cohen's *d* |
| 7 | Pearson correlation: theta ~ alpha |
| 8 | Visualisations (3 figures) |
| 9 | Final summary table and CSV export |

---

## 1. Data Loading

Two input files are loaded:

- **`results_all_subjects.csv`** — one row per subject, produced by the preprocessing pipeline. Contains band power in dB for each of the five canonical bands (δ, θ, α, β, γ).
- **`participants.tsv`** — the BIDS-mandated metadata file at the dataset root. Provides the `group` column (`expert` / `novice`), which is the ground truth for group assignment. Subject IDs are normalised to zero-padded three-digit strings (e.g. `'1'` → `'001'`) to ensure consistent string matching across both files.

---

## 2. Merge, Group Assignment, and Data Cleaning

The two DataFrames are joined on the shared `subject` key using an **inner join**:

$$\text{df} = \text{df\_results} \bowtie_{\text{subject}} \text{df\_participants}[\texttt{subject, group, age, gender}]$$

An inner join retains only subjects present in both files, automatically handling edge cases such as preprocessing failures or subjects with metadata but no EEG data. The `group` label is sourced directly from `participants.tsv` and is never hardcoded.

After merging, **sub-007** is excluded (anomalous raw PSD identified during visual QC in the preprocessing notebook), leaving a final analytic sample of **N = 23** (11 experts, 12 novices).

---

## 3. Descriptive Statistics by Group

For each group and each frequency band, the following statistics are computed using pandas `.describe()` with `ddof=1`:

| Statistic | Formula | Meaning |
|-----------|---------|---------|
| Mean | $\bar{x} = \frac{1}{n}\sum_{i=1}^n x_i$ | Central tendency |
| Standard deviation | $s = \sqrt{\frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2}$ | Spread (unbiased, ddof=1) |
| Median | $x_{(\lceil n/2 \rceil)}$ | Robust centre |
| Min / Max | — | Observed range |

The unbiased estimator (dividing by $n-1$) is used throughout because the groups are samples, not populations.

---

## 4. Welch's t-test

### Rationale

Student's t-test assumes homoscedasticity (equal population variances). With small and unequal group sizes ($n_E = 11$, $n_N = 12$) and visibly different standard deviations between groups (e.g. alpha: $s_E \approx 2.9$ dB vs $s_N \approx 5.4$ dB), this assumption cannot be taken for granted. Welch's test relaxes it by using separate variance estimates for each group.

### Test statistic

$$t = \frac{\bar{x}_E - \bar{x}_N}{\sqrt{\dfrac{s_E^2}{n_E} + \dfrac{s_N^2}{n_N}}}$$

### Degrees of freedom — Welch–Satterthwaite approximation

Because the two variances are estimated separately, the degrees of freedom are not simply $n_E + n_N - 2$. The Welch–Satterthwaite equation gives a non-integer approximation:

$$\nu = \frac{\left(\dfrac{s_E^2}{n_E} + \dfrac{s_N^2}{n_N}\right)^2}{\dfrac{(s_E^2/n_E)^2}{n_E-1} + \dfrac{(s_N^2/n_N)^2}{n_N-1}}$$

The test is performed with `scipy.stats.ttest_ind(equal_var=False)`, which implements both the statistic and the Satterthwaite df internally. Raw p-values are collected for FDR correction in the next step.

---

## 5. Benjamini–Hochberg FDR Correction

### The multiple comparisons problem

Five independent t-tests (one per frequency band) are performed. Without correction, the probability of at least one false positive by chance is:

$$P(\geq 1 \text{ false positive}) = 1 - (1-\alpha)^m = 1 - 0.95^5 \approx 0.226$$

### The Benjamini–Hochberg (BH) procedure

Rather than controlling the familywise error rate (FWER, as Bonferroni does), BH controls the **False Discovery Rate** — the expected proportion of false positives among all rejected hypotheses. The procedure is more powerful than Bonferroni and is the standard choice when $m$ is small and tests are independent or positively correlated.

**Algorithm** (implemented from scratch with NumPy, no external packages):

1. Sort the $m$ raw p-values in ascending order: $p_{(1)} \leq p_{(2)} \leq \cdots \leq p_{(m)}$
2. For each rank $k$, compute the BH critical value: $p_{(k)} \leq \dfrac{k}{m} \alpha$
3. Find the largest $k^*$ such that $p_{(k^*)} \leq \dfrac{k^*}{m}\alpha$; reject all hypotheses with rank $\leq k^*$
4. Compute adjusted p-values using `np.minimum.accumulate` on the reversed sorted array to ensure monotonicity

### Results

After BH-FDR correction at $\alpha = 0.05$, **no frequency band survived** — all adjusted p-values equal 0.660. This is consistent with the small sample size, high within-group variability, and small observed effect sizes.

---

## 6. Effect Size: Cohen's *d*

### Rationale

The p-value answers *"Is this difference unlikely under $H_0$?"* but says nothing about magnitude. With small samples ($n \approx 11$), even a large true effect can yield a non-significant p-value due to low statistical power. Cohen's $d$ quantifies the group difference in standardised, sample-size-independent units.

### Formula

$$d = \frac{\bar{x}_E - \bar{x}_N}{s_{\text{pooled}}}, \qquad s_{\text{pooled}} = \sqrt{\frac{(n_E-1)\,s_E^2 + (n_N-1)\,s_N^2}{n_E + n_N - 2}}$$

The pooled SD weights each group's variance by its degrees of freedom. Cohen's conventional benchmarks are applied: $|d| < 0.2$ = negligible, $0.2$–$0.5$ = small, $0.5$–$0.8$ = medium, $> 0.8$ = large.

### Results

| Band | Cohen's *d* | Interpretation |
|------|------------|----------------|
| Delta | 0.237 | Small |
| Theta | −0.211 | Small |
| Alpha | 0.336 | Small |
| Beta | −0.187 | Small |
| Gamma | −0.397 | Small |

All effect sizes are small ($|d| < 0.4$), consistent with the null FDR results. The direction of the gamma effect (novices > experts) is noteworthy but underpowered.

---

## 7. Pearson Correlation: Theta ~ Alpha

### Rationale

In the absence of subjective meditation-depth ratings, the cross-subject linear relationship between theta-band and alpha-band power is examined as a proxy measure of shared neural dynamics. Both bands are theoretically linked to meditative states: theta to hippocampal–frontal memory encoding, alpha to cortical attentional gating.

> ⚠️ This is a **between-subjects** correlation: each data point is one participant's mean band power over the entire session. It does not capture within-subject dynamics.

### Inference

$$r = \frac{\sum_i (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_i(x_i-\bar{x})^2}\sqrt{\sum_i(y_i-\bar{y})^2}}$$

A 95% confidence interval is computed via Fisher's z-transformation:

$$z_r = \text{arctanh}(r), \qquad \text{SE}(z_r) = \frac{1}{\sqrt{n-3}}, \qquad \text{CI} = \left[\tanh\!\left(z_r \mp 1.96 \cdot \text{SE}\right)\right]$$

### Result

$$r = 0.763, \quad p < 0.001, \quad 95\%\,\text{CI} = [0.52,\, 0.89]$$

A strong positive correlation between theta and alpha power is observed across all subjects ($n = 23$), suggesting that participants with overall higher spectral power in the theta range also show proportionally higher alpha power, regardless of group membership.

---

## 8. Visualisations

Three figures are produced and saved as PNG files:

### 8a — Boxplot with individual data points (`group_comparison_boxplot.png`)

One panel per frequency band (5 panels). Each panel shows side-by-side boxplots for expert and novice groups. Individual subject values are overlaid as jittered scatter points to expose the full distribution rather than just summary statistics. A significance bracket is drawn above bands where $p^{\text{BH}} < 0.05$ (none in this dataset).

### 8b — Theta ~ Alpha scatter plot with OLS regression line (`theta_alpha_correlation.png`)

Each point represents one subject, coloured by group (blue = expert, red = novice). Each point is labelled with the subject ID. An ordinary least-squares regression line is fitted across all subjects. The Pearson correlation, 95% CI, and p-value computed in Section 7 are annotated in the corner of the plot.

### 8c — Cohen's *d* bar chart (`cohens_d_barplot.png`)

One bar per frequency band. Bar colour encodes direction: blue if experts have higher power ($d > 0$), red if novices do ($d < 0$). Dashed horizontal reference lines mark Cohen's small (0.2), medium (0.5), and large (0.8) thresholds on both sides of zero.

---

## 9. Summary Table and Export

All results are assembled into a single publication-ready table and exported to `statistical_results.csv`. Columns include:

| Column | Content |
|--------|---------|
| Band | Frequency band name and range |
| M Expert / M Novice (dB) | Group means |
| SD Expert / SD Novice | Group standard deviations |
| *t* | Welch t-statistic |
| df (Welch–Satterthwaite) | Non-integer degrees of freedom |
| p (raw) | Uncorrected two-tailed p-value |
| p (BH-FDR) | BH-adjusted p-value |
| Cohen's *d* | Standardised effect size |
| Effect size | Verbal label (negligible / small / medium / large) |
| Sig. | Significance marker (`*` / `n.s.`) |

### Final Results

| Band | M Expert | M Novice | *t* | df | p (BH-FDR) | Cohen's *d* | Sig. |
|------|----------|----------|-----|----|------------|-------------|------|
| Delta | −114.17 | −115.05 | 0.553 | 15.2 | 0.660 | 0.237 | n.s. |
| Theta | −119.12 | −118.29 | −0.512 | 20.3 | 0.660 | −0.211 | n.s. |
| Alpha | −114.84 | −116.30 | 0.826 | 17.2 | 0.660 | 0.336 | n.s. |
| Beta | −124.43 | −123.60 | −0.446 | 20.6 | 0.660 | −0.187 | n.s. |
| Gamma | −130.34 | −128.80 | −0.974 | 17.8 | 0.660 | −0.397 | n.s. |

No significant group differences were observed after FDR correction. The most plausible explanations are the small analytic sample ($N = 23$, ~11 per group), high within-group variability in band power (SD 2–5 dB), and the presence of sub-003 as a high-power outlier within the expert group.
