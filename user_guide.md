# Power Panel — User Guide

## 1. Overview

Power Panel is a web application for conducting statistical power analyses for commonly used panel data designs. It implements the methods developed in:

> Schochet, P.Z. (2022). *Statistical Power for Estimating Treatment Effects Using Difference-in-Differences and Comparative Interrupted Time Series Estimators with Variation in Treatment Timing.* Journal of Educational and Behavioral Statistics, 48(6), 713–751.

The app supports three research designs, two analysis modes, and a scenario comparison feature for exploring how power varies across parameter assumptions.

---

## 2. Getting Started

### Installation

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

### Layout

The interface has two main areas:

- **Sidebar** (left): All input parameters, organized into five collapsible sections.
- **Main area** (center/right): Three tabs — *Results*, *Scenario Comparison*, and *Help*.

---

## 3. Sidebar Parameters

### 3.1 Design Configuration

#### Design Type

| Design | Description | When to use |
|--------|-------------|-------------|
| **DID** | Difference-in-Differences | Treatment and comparison groups observed before and after treatment. The simplest design — compares changes over time between groups. |
| **CITS** | Comparative Interrupted Time Series | Like DID but also models time trends (slopes). Requires more time periods but can detect trend changes. Needs a comparison group. |
| **ITS** | Interrupted Time Series | Models the time trend for treated units only — no comparison group required. Useful when no valid comparison group exists. |

#### Linear Model Specification (CITS / ITS only)

These options control how the time trend is modeled:

- **Fully interacted (different slopes)**: Allows separate pre- and post-treatment slopes for treatment and comparison groups. Most flexible but requires the most data.
- **Common slopes**: Assumes treatment and comparison groups share the same time trend slope. More efficient if the assumption holds.
- **Discrete post-period indicators**: Replaces the post-period linear trend with separate indicator variables for each post-treatment time period (similar to a DID-style specification). Useful when effects may not follow a linear trend.

#### Analysis Mode

- **Calculate MDE given sample size**: You provide the number of clusters; the app computes the minimum detectable effect size (MDE) at each power level.
- **Calculate required clusters given target MDE**: You provide a target MDE; the app computes how many clusters are needed at each power level.

#### Panel Data Type

- **Cross-sectional**: Different individuals are observed in each time period (repeated cross-sections).
- **Longitudinal**: The same individuals are tracked across all time periods (panel data).

The choice affects how individual-level autocorrelation is handled. In longitudinal panels, within-person correlation (phi) reduces residual variance and improves power.

---

### 3.2 Time Periods & Timing Groups

#### Number of Time Periods

The total number of observation periods (pre + post treatment). More time periods generally improve power, especially for CITS/ITS designs that model trends.

**Minimum requirements:**
- DID: at least 2 periods
- CITS/ITS with slopes (fully interacted or common slopes): at least 6 periods
- CITS/ITS discrete: at least 4 periods

#### Evenly Spaced Intervals

- **On** (default): Time periods are equally spaced (e.g., years 1, 2, 3, …).
- **Off**: Enter custom intervals (e.g., `1, 2, 3, 5, 8` for uneven spacing). This matters for designs that model time trends, since the spacing affects how well slopes are estimated.

#### Number of Timing Groups

A timing group is a set of clusters that begin treatment at the same time. This is how the app handles **staggered treatment adoption** — when different clusters start treatment in different periods.

- **1 group**: All treated clusters start at the same time (simple design).
- **2+ groups**: Different clusters start treatment at different times (staggered rollout).

#### Configuring Each Timing Group

For each timing group, you specify:

- **Start period**: The first period in which the group receives treatment. Must be between 2 and the total number of periods (period 1 is always pre-treatment). Start periods must be in increasing order across groups.

Depending on the analysis mode:

- **MDE mode**: Enter the number of **treatment clusters** (and **comparison clusters** for DID/CITS) assigned to this group.
- **Required clusters mode**: Enter the **treatment share** (and **comparison share** for DID/CITS) — the proportion of total clusters allocated to this group. Shares must sum to 1.0 across all groups.

#### Post-Treatment Period for Analysis

Choose what post-treatment outcome to analyze:

- **Average across all post-periods**: The estimated effect is the average treatment effect across all post-treatment time points. Generally provides the most statistical power.
- **Single post-period time point**: The estimated effect at a specific time. Choose between:
  - *At a specific time period*: e.g., "the effect at period 8"
  - *At a specific treatment exposure point*: e.g., "the effect after 3 periods of exposure" (adjusts for staggered timing — each group's effect is measured at the same exposure duration)

---

### 3.3 Sample Sizes

#### Individuals per Cluster per Period

The number of individuals observed in each cluster at each time period. Larger clusters reduce individual-level sampling variance but provide diminishing returns once cluster-level variance dominates (high ICC).

#### Mode-specific inputs

- **MDE mode**: The cluster counts you entered per timing group are summed and displayed as totals.
- **Required clusters mode**: Enter the **proportion of all clusters that are treatment** (e.g., 0.50 for an equal split). This proportion applies to the total cluster count that the app solves for.

---

### 3.4 Error Structure

#### Autocorrelation Structure

Controls how outcomes are correlated within clusters over time:

- **AR(1)**: First-order autoregressive — correlation decays geometrically with time distance. A correlation at lag 1 of rho implies a correlation of rho^k at lag k. Common in many social science settings.
- **Constant**: Equal correlation (rho) between all time periods regardless of distance. Equivalent to a compound symmetry structure.
- **None**: No within-cluster autocorrelation over time.

#### Cluster Autocorrelation (rho)

The correlation of cluster-level means across time periods. Range: −0.99 to 0.99. Higher values (e.g., 0.40–0.80) are common when clusters are stable entities like schools or clinics. Positive autocorrelation *increases* variance in DID-type designs because the "difference" does not remove as much noise.

#### Individual Autocorrelation (phi) — Longitudinal panels only

The within-person correlation across time periods. Only relevant for longitudinal data where the same individuals are tracked. Higher values mean person-level outcomes are more predictable over time, which *reduces* residual variance and improves power.

#### Intraclass Correlation (ICC)

The proportion of total variance that is between clusters (as opposed to within clusters). Range: 0 to 1.

- Low ICC (0.01–0.05): Most variation is within clusters. Common for large clusters (e.g., large schools).
- High ICC (0.10–0.30): Substantial between-cluster variation. Common for small clusters or when clusters differ substantially.

Higher ICC increases the required sample size because cluster-level sampling variance dominates.

---

### 3.5 Precision & Testing

#### R-squared from Covariates (R²yx)

The proportion of outcome variance explained by model covariates (e.g., baseline demographics, pre-treatment measures). Higher values reduce residual variance and improve power. Set to 0 if no covariates are included.

#### Treatment-Covariate Correlation (R²tx)

The proportion of treatment status variance explained by covariates. If covariates predict treatment assignment, this *reduces* the effective variation in treatment, which *hurts* power. Usually small in randomized designs. Set to 0 for randomized experiments.

#### Design Effect (Weighting)

A multiplicative adjustment for survey weights or complex sampling. Default is 1.0 (no adjustment). Values > 1 increase variance and reduce power.

#### Significance Level (alpha)

The probability of a Type I error (false positive). Standard value is 0.05.

#### Two-Tailed Test

- **On** (default): Tests for effects in either direction. Use for most applications.
- **Off**: One-tailed test — only tests for effects in one direction. More powerful but only appropriate when the direction of the effect is known a priori.

#### Power Range

The range of statistical power levels for which to compute results (MDE or required clusters). Default: 0.60 to 0.90. Power is the probability of detecting a true effect of the given size.

#### Power Step

The increment between power levels in the output table (0.01, 0.05, or 0.10).

---

## 4. Results Tab

When all inputs are valid, the Results tab displays:

1. **Design label**: A summary of the selected design (e.g., "DID — Cross-Sectional — AR(1)").

2. **Results table**: A table with one row per power level. Columns are:
   - *Power*: The statistical power level.
   - *MDE* (mode 1): The minimum detectable effect size in standard deviation units. Smaller is better — it means you can detect smaller effects.
   - *Required Clusters* (mode 2): The total number of clusters needed. This is the sum across treatment and comparison groups.

3. **Line chart**: An interactive Plotly chart showing the same data. Hover over points for exact values. You can zoom, pan, and download the chart as an image.

4. **Input Parameters Summary** (collapsible): A text summary of all input parameters used for the computation, for documentation purposes.

5. **Download CSV**: Export the results table as a CSV file.

### Interpreting MDE Values

The MDE is expressed in standard deviation units (Cohen's d). Common benchmarks:

| MDE | Interpretation |
|-----|---------------|
| 0.10 | Very small effect — requires very large samples |
| 0.20 | Small effect — typical target for well-powered studies |
| 0.30 | Small-to-medium effect |
| 0.50 | Medium effect |

Lower MDE values indicate greater statistical precision but require more clusters or individuals.

---

## 5. Scenario Comparison Tab

This feature lets you explore how MDE (or required clusters) varies across different assumptions about ICC and R-squared values. It is especially useful for:

- **Grant proposals**: Show reviewers that the study is well-powered under a range of plausible assumptions.
- **Sensitivity analysis**: Understand which parameters have the largest impact on power.
- **Study planning**: Identify the parameter combinations where your design becomes underpowered.

### How to use it

1. **ICC values**: Enter a comma-separated list of ICC values to compare (e.g., `0.01, 0.05, 0.10, 0.15, 0.20`).
2. **R²yx values**: Enter a comma-separated list of R-squared values (e.g., `0.00, 0.10, 0.20, 0.30, 0.50`).
3. **Fixed power level**: Choose a single power level for the comparison (e.g., 0.80).
4. Click **Generate Scenario Table**.

### Output

A table with:
- **Rows**: One per ICC value.
- **Columns**: One per R²yx value.
- **Cells**: The MDE (or required clusters) at the specified power level, using all other parameters from the sidebar.

The table is designed for easy **copy-paste** into Word, Excel, or other documents. You can also download it as CSV.

### Example interpretation

Suppose the table shows MDE = 0.18 at ICC = 0.05 and R²yx = 0.20, but MDE = 0.31 at ICC = 0.15 and R²yx = 0.00. This tells you that:
- With moderate covariates and low ICC, you can detect effects as small as 0.18 SD.
- Without covariates and with higher ICC, the smallest detectable effect rises to 0.31 SD.
- Adding covariates (increasing R²yx) has a substantial benefit for power.

---

## 6. Key Concepts

### What is a Minimum Detectable Effect Size (MDE)?

The MDE is the smallest true treatment effect that your study design can reliably detect at a given power level and significance level. It depends on:

- Sample size (clusters and individuals)
- Number of time periods
- ICC and autocorrelation
- Covariate precision gains
- Design type

The formula is:

```
MDE = Factor(alpha, power, df) × sqrt(Var(beta_hat))
```

where `Factor` combines the critical values from the t-distribution, and `Var(beta_hat)` is the sampling variance of the treatment effect estimator.

### Staggered Treatment Timing

In many real-world evaluations, not all units start treatment at the same time. The app handles this through **timing groups**: each group has its own treatment start period, and the variance formulas account for the different pre/post splits across groups.

More timing groups with different start times generally *reduce* power because each group contributes fewer effective observations to the trend estimation. However, staggered designs are often necessary for practical or ethical reasons.

### Cross-Sectional vs. Longitudinal

- **Cross-sectional panels**: New individuals in each period. Within-person correlation is not relevant. Power depends mainly on cluster-level autocorrelation (rho) and ICC.
- **Longitudinal panels**: Same individuals across periods. Within-person correlation (phi) adds information that reduces variance. Generally provides better power than cross-sectional designs for the same sample sizes.

---

## 7. Tips for Practitioners

1. **Start with the scenario comparison.** Before committing to a sample size, generate a grid to see how sensitive your MDE is to ICC and R²yx assumptions. These are often the hardest parameters to pin down in advance.

2. **ICC matters more than you think.** A jump from ICC = 0.05 to ICC = 0.15 can double your required sample size. If possible, use pilot data or prior studies to estimate ICC.

3. **Covariates are free power.** Adding baseline covariates (increasing R²yx) reduces MDE without adding clusters. Even modest R²yx = 0.20 can meaningfully improve power.

4. **More time periods help CITS/ITS more than DID.** For DID, additional periods help mainly through autocorrelation. For CITS/ITS, more periods directly improve slope estimation.

5. **Consider the post-period choice.** Averaging across all post-periods is usually more powerful than looking at a single point. But if you care about the effect at a specific time (e.g., 2 years post-treatment), use the specific period option.

6. **Longitudinal > Cross-sectional for power** (all else equal), because within-person tracking adds information. But longitudinal panels are harder and more expensive to collect.

7. **Two-tailed vs. one-tailed.** Use two-tailed unless you have a strong a priori directional hypothesis. One-tailed tests are more powerful but can miss effects in the unexpected direction.

---

## 8. Troubleshooting

| Problem | Solution |
|---------|----------|
| "No results computed" | Check that all inputs are valid. Look for red error messages above the results. |
| Very large MDE values | Your sample may be too small. Try increasing clusters, individuals per cluster, or R²yx. |
| Required clusters is extremely high | The target MDE may be too ambitious for your design. Consider a larger effect size or adding covariates. |
| App is slow | The scenario comparison runs many computations. Reduce the number of ICC or R² values in the grid. |
| Validation errors about start times | Start times must be in increasing order, between 2 and the total number of periods. |

---

## 9. References

- Schochet, P.Z. (2022). Statistical Power for Estimating Treatment Effects Using Difference-in-Differences and Comparative Interrupted Time Series Estimators with Variation in Treatment Timing. *Journal of Educational and Behavioral Statistics*, 48(6), 713–751. https://doi.org/10.3102/10769986221112426

---

*Power Panel is open source under the MIT license. Source code: [github.com/facusabino/power-analysis-app](https://github.com/facusabino/power-analysis-app)*
