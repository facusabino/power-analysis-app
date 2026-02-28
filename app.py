"""
Power Panel: Streamlit app for power analysis of panel data designs.
Based on Schochet (JEBS, 2022).
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

from config import (
    PowerInputs, DesignType, CITSSpec, AnalysisMode,
    PanelType, AutocorrStructure, PostPeriodType, SpecificPostUnit,
)
from validation import validate_inputs
from power_engine import run_power_analysis
from scenarios import generate_scenario_grid

# --- Page config ---
st.set_page_config(page_title="Power Panel", layout="wide")
st.title("Power Analysis for Panel Data Designs")
st.caption("Based on Schochet (Journal of Educational and Behavioral Statistics, 2022)")

# =====================================================================
# SIDEBAR
# =====================================================================

with st.sidebar:
    st.header("Parameters")

    # --- Section 1: Design Configuration ---
    with st.expander("Design Configuration", expanded=True):
        design_choice = st.radio(
            "Design type",
            ["Difference-in-Differences (DID)",
             "Comparative Interrupted Time Series (CITS)",
             "Interrupted Time Series (ITS)"],
        )
        design_map = {
            "Difference-in-Differences (DID)": DesignType.DID,
            "Comparative Interrupted Time Series (CITS)": DesignType.CITS,
            "Interrupted Time Series (ITS)": DesignType.ITS,
        }
        design_type = design_map[design_choice]

        cits_spec = CITSSpec.FULLY_INTERACTED
        if design_type in (DesignType.CITS, DesignType.ITS):
            spec_choice = st.radio(
                "Linear model specification",
                ["Fully interacted (different slopes)",
                 "Common slopes",
                 "Discrete post-period indicators"],
            )
            spec_map = {
                "Fully interacted (different slopes)": CITSSpec.FULLY_INTERACTED,
                "Common slopes": CITSSpec.COMMON_SLOPES,
                "Discrete post-period indicators": CITSSpec.DISCRETE,
            }
            cits_spec = spec_map[spec_choice]

        mode_choice = st.radio(
            "Analysis mode",
            ["Calculate MDE given sample size",
             "Calculate required clusters given target MDE"],
        )
        analysis_mode = AnalysisMode.CALC_MDE if "MDE given" in mode_choice else AnalysisMode.CALC_CLUSTERS

        target_mde = 0.20
        if analysis_mode == AnalysisMode.CALC_CLUSTERS:
            target_mde = st.number_input("Target MDE", value=0.20, min_value=0.01, step=0.01)

        panel_choice = st.radio(
            "Panel data type",
            ["Cross-sectional (different people over time)",
             "Longitudinal (same people over time)"],
        )
        panel_type = PanelType.CROSS_SECTIONAL if "Cross" in panel_choice else PanelType.LONGITUDINAL

    # --- Section 2: Time Periods & Timing Groups ---
    with st.expander("Time Periods & Timing Groups"):
        n_time_periods = st.number_input("Number of time periods", min_value=2, value=10, step=1)

        evenly_spaced = st.toggle("Evenly spaced intervals", value=True)
        time_intervals = None
        if not evenly_spaced:
            ti_str = st.text_input(
                "Time intervals (comma-separated integers, earliest to most recent)",
                value=",".join(str(i) for i in range(1, n_time_periods + 1)),
            )
            try:
                time_intervals = [int(x.strip()) for x in ti_str.split(",")]
            except ValueError:
                st.error("Time intervals must be comma-separated integers.")
                time_intervals = list(range(1, n_time_periods + 1))

        n_timing_groups = st.number_input("Number of timing groups", min_value=1, value=2, step=1)

        st.markdown("**Configure each timing group:**")

        start_times = []
        mtk_values = []
        mck_values = []
        rtk_values = []
        rck_values = []

        n_cols = min(n_timing_groups, 3)
        cols = st.columns(n_cols)

        for k in range(n_timing_groups):
            with cols[k % n_cols]:
                st.markdown(f"**Group {k+1}**")
                default_start = min(4 + k * 2, n_time_periods)
                st_val = st.number_input(
                    f"Start period",
                    min_value=2, max_value=n_time_periods,
                    value=default_start,
                    key=f"start_{k}",
                )
                start_times.append(st_val)

                if analysis_mode == AnalysisMode.CALC_MDE:
                    mt_val = st.number_input(
                        f"Treatment clusters",
                        min_value=2, value=10, key=f"mt_{k}",
                    )
                    mtk_values.append(mt_val)

                    if design_type != DesignType.ITS:
                        mc_val = st.number_input(
                            f"Comparison clusters",
                            min_value=2, value=10, key=f"mc_{k}",
                        )
                        mck_values.append(mc_val)
                else:
                    default_share = round(1.0 / n_timing_groups, 2)
                    rtk_val = st.number_input(
                        f"Treatment share",
                        min_value=0.01, max_value=1.0,
                        value=default_share, step=0.01,
                        key=f"rtk_{k}",
                    )
                    rtk_values.append(rtk_val)

                    if design_type != DesignType.ITS:
                        rck_val = st.number_input(
                            f"Comparison share",
                            min_value=0.01, max_value=1.0,
                            value=default_share, step=0.01,
                            key=f"rck_{k}",
                        )
                        rck_values.append(rck_val)

        # Post-period type
        post_choice = st.radio(
            "Post-treatment period for analysis",
            ["Average across all post-periods", "Single post-period time point"],
        )
        post_period_type = PostPeriodType.AVERAGE if "Average" in post_choice else PostPeriodType.SPECIFIC

        specific_post_unit = SpecificPostUnit.TIME_PERIOD
        q_time = 6
        l_time = 6
        if post_period_type == PostPeriodType.SPECIFIC:
            unit_choice = st.radio(
                "Measurement unit",
                ["At a specific time period", "At a specific treatment exposure point"],
            )
            specific_post_unit = SpecificPostUnit.TIME_PERIOD if "time period" in unit_choice else SpecificPostUnit.EXPOSURE_POINT

            if specific_post_unit == SpecificPostUnit.TIME_PERIOD:
                q_time = st.number_input("Specific time point", min_value=1, value=6, step=1)
            else:
                l_time = st.number_input("Specific exposure point", min_value=1, value=6, step=1)

    # --- Section 3: Sample Sizes ---
    with st.expander("Sample Sizes"):
        n_per_cluster = st.number_input("Individuals per cluster per period", min_value=1, value=100, step=1)

        if analysis_mode == AnalysisMode.CALC_MDE:
            total_mt = sum(mtk_values) if mtk_values else 20
            st.metric("Total treatment clusters", total_mt)
            if design_type != DesignType.ITS:
                total_mc = sum(mck_values) if mck_values else 20
                st.metric("Total comparison clusters", total_mc)
        else:
            rt = 0.50
            if design_type != DesignType.ITS:
                rt = st.number_input("Proportion of all clusters that are treatment",
                                     min_value=0.01, max_value=0.99, value=0.50, step=0.01)

    # --- Section 4: Error Structure ---
    with st.expander("Error Structure"):
        autocorr_choice = st.radio("Autocorrelation structure", ["AR1", "Constant", "None"], horizontal=True)
        autocorr_map = {"AR1": AutocorrStructure.AR1, "Constant": AutocorrStructure.CONSTANT, "None": AutocorrStructure.NONE}
        autocorr_structure = autocorr_map[autocorr_choice]

        rho = 0.0
        phi = 0.0
        if autocorr_structure != AutocorrStructure.NONE:
            rho = st.slider("Cluster autocorrelation (rho)", -0.99, 0.99, 0.40, 0.01)
            if panel_type == PanelType.LONGITUDINAL:
                phi = st.slider("Individual autocorrelation (phi)", -0.99, 0.99, 0.40, 0.01)

        icc = st.slider("Intraclass correlation (ICC)", 0.0, 1.0, 0.05, 0.01)

    # --- Section 5: Precision & Testing ---
    with st.expander("Precision & Testing"):
        r2yx = st.slider("R-squared from covariates (R2yx)", 0.0, 0.99, 0.0, 0.01)
        r2tx = st.slider("Treatment-covariate correlation (R2tx)", 0.0, 0.99, 0.0, 0.01)
        deff_wgt = st.number_input("Design effect (weighting)", min_value=1.0, value=1.0, step=0.1)
        alpha = st.number_input("Significance level (alpha)", min_value=0.01, max_value=0.50, value=0.05, step=0.01)
        two_tailed = st.toggle("Two-tailed test", value=True)
        power_range = st.slider("Power range", 0.50, 0.99, (0.60, 0.90), 0.01)
        power_step = st.select_slider("Power step", options=[0.01, 0.05, 0.10], value=0.05)


# =====================================================================
# BUILD INPUTS
# =====================================================================

mt = sum(mtk_values) if mtk_values else 20
mc = sum(mck_values) if mck_values else 20

inputs = PowerInputs(
    design_type=design_type,
    cits_spec=cits_spec,
    analysis_mode=analysis_mode,
    panel_type=panel_type,
    n_time_periods=n_time_periods,
    evenly_spaced=evenly_spaced,
    time_intervals=time_intervals,
    n_timing_groups=n_timing_groups,
    start_times=start_times,
    post_period_type=post_period_type,
    specific_post_unit=specific_post_unit,
    q_time=q_time,
    l_time=l_time,
    mt=mt,
    mtk=mtk_values if mtk_values else [10, 10],
    mc=mc,
    mck=mck_values if mck_values else [10, 10],
    target_mde=target_mde,
    rt=rt if analysis_mode == AnalysisMode.CALC_CLUSTERS else 0.50,
    rtk=rtk_values if rtk_values else [0.50, 0.50],
    rck=rck_values if rck_values else [0.50, 0.50],
    n_per_cluster=n_per_cluster,
    autocorr_structure=autocorr_structure,
    rho=rho,
    phi=phi,
    icc=icc,
    r2yx=r2yx,
    r2tx=r2tx,
    deff_wgt=deff_wgt,
    alpha=alpha,
    two_tailed=two_tailed,
    power_min=power_range[0],
    power_max=power_range[1],
    power_step=power_step,
)


# =====================================================================
# MAIN AREA WITH TABS
# =====================================================================

tab_results, tab_scenarios, tab_help = st.tabs(["Results", "Scenario Comparison", "Help"])

# --- Results Tab ---
with tab_results:
    validation = validate_inputs(inputs)

    if not validation.is_valid:
        for err in validation.errors:
            st.error(err)
    else:
        result = run_power_analysis(inputs)

        if 'error' in result:
            st.error(result['error'])
        elif result['results']:
            st.subheader(result['design_label'])

            df = pd.DataFrame(result['results'])

            if analysis_mode == AnalysisMode.CALC_MDE:
                df.columns = ["Power", "MDE"]
                fmt = {"Power": "{:.2f}", "MDE": "{:.4f}"}
                y_col = "MDE"
                y_label = "Minimum Detectable Effect Size"
            else:
                df.columns = ["Power", "Required Clusters"]
                fmt = {"Power": "{:.2f}", "Required Clusters": "{:,.0f}"}
                y_col = "Required Clusters"
                y_label = "Required Number of Clusters"

            col_table, col_chart = st.columns([1, 1.5])

            with col_table:
                st.dataframe(
                    df.style.format(fmt),
                    use_container_width=False,
                    hide_index=True,
                    height=min(len(df) * 40 + 40, 600),
                )

            with col_chart:
                fig = px.line(
                    df, x="Power", y=y_col,
                    markers=True,
                    labels={"Power": "Power Level", y_col: y_label},
                )
                fig.update_traces(line_color="#00aff5", line_width=2.5)
                fig.update_layout(
                    margin=dict(l=20, r=20, t=30, b=20),
                    height=400,
                )
                st.plotly_chart(fig, use_container_width=True)

            # Input summary
            with st.expander("Input Parameters Summary"):
                panel_label = "Cross-Sectional" if panel_type == PanelType.CROSS_SECTIONAL else "Longitudinal"
                summary = f"""
- **Design**: {result['design_label']}
- **Time periods**: {n_time_periods}, **Timing groups**: {n_timing_groups}
- **Start times**: {start_times}
- **Panel type**: {panel_label}
- **Autocorrelation**: {autocorr_choice}, rho={rho:.2f}"""
                if panel_type == PanelType.LONGITUDINAL:
                    summary += f", phi={phi:.2f}"
                summary += f"""
- **ICC**: {icc:.2f}
- **R2yx**: {r2yx:.2f}, **R2tx**: {r2tx:.2f}
- **Design effect**: {deff_wgt:.1f}
- **Alpha**: {alpha:.2f}, **{'Two' if two_tailed else 'One'}-tailed test**
- **Individuals per cluster**: {n_per_cluster}"""
                if analysis_mode == AnalysisMode.CALC_MDE:
                    summary += f"\n- **Treatment clusters**: {mt}, **Comparison clusters**: {mc}"
                else:
                    summary += f"\n- **Target MDE**: {target_mde:.2f}, **Treatment proportion**: {rt:.2f}"
                st.markdown(summary)

            # CSV download
            csv = df.to_csv(index=False)
            st.download_button("Download results as CSV", csv, "power_results.csv", "text/csv")

        else:
            st.warning("No results computed. Check your inputs.")


# --- Scenario Comparison Tab ---
with tab_scenarios:
    st.subheader("Scenario Comparison")
    st.markdown("Compare MDE across different ICC and R-squared values. Configure the grid below.")

    col1, col2, col3 = st.columns(3)
    with col1:
        icc_str = st.text_input(
            "ICC values (comma-separated)",
            value="0.01, 0.05, 0.10, 0.15, 0.20",
            key="scenario_icc",
        )
    with col2:
        r2_str = st.text_input(
            "R2yx values (comma-separated)",
            value="0.00, 0.10, 0.20, 0.30, 0.50",
            key="scenario_r2",
        )
    with col3:
        scenario_power = st.number_input(
            "Fixed power level",
            value=0.80, min_value=0.50, max_value=0.99, step=0.01,
            key="scenario_power",
        )

    run_scenarios = st.button("Generate Scenario Table")

    if run_scenarios:
        try:
            icc_list = [float(x.strip()) for x in icc_str.split(",")]
            r2_list = [float(x.strip()) for x in r2_str.split(",")]
        except ValueError:
            st.error("Please enter valid comma-separated numbers.")
            icc_list = []
            r2_list = []

        if icc_list and r2_list:
            validation = validate_inputs(inputs)
            if not validation.is_valid:
                for err in validation.errors:
                    st.error(err)
            else:
                with st.spinner("Computing scenarios..."):
                    scenario_df = generate_scenario_grid(inputs, icc_list, r2_list, scenario_power)

                if analysis_mode == AnalysisMode.CALC_MDE:
                    st.markdown(f"**MDE values at power = {scenario_power:.2f}**")
                    st.dataframe(
                        scenario_df.style.format("{:.4f}"),
                        use_container_width=True,
                    )
                else:
                    st.markdown(f"**Required clusters at power = {scenario_power:.2f}**")
                    st.dataframe(
                        scenario_df.style.format("{:,.0f}"),
                        use_container_width=True,
                    )

                csv = scenario_df.to_csv()
                st.download_button(
                    "Download scenario table as CSV",
                    csv,
                    "scenario_comparison.csv",
                    "text/csv",
                )


# --- Help Tab ---
with tab_help:
    st.subheader("User Guide")

    # Download link for the full guide
    import pathlib
    guide_path = pathlib.Path(__file__).parent / "user_guide.md"
    if guide_path.exists():
        guide_text = guide_path.read_text(encoding="utf-8")
        st.download_button(
            "Download full user guide (Markdown)",
            guide_text,
            "Power_Panel_User_Guide.md",
            "text/markdown",
        )

    st.markdown("---")

    with st.expander("Design Types", expanded=True):
        st.markdown("""
| Design | Description |
|--------|-------------|
| **DID** | Difference-in-Differences — compares changes over time between treatment and comparison groups. Simplest design. |
| **CITS** | Comparative Interrupted Time Series — also models time trends (slopes). Requires more time periods but can detect trend changes. Needs a comparison group. |
| **ITS** | Interrupted Time Series — models the time trend for treated units only. No comparison group required. |

**CITS/ITS model specifications:**
- **Fully interacted**: Separate pre/post slopes for each group. Most flexible, most data-hungry.
- **Common slopes**: Shared time trend across groups. More efficient if the assumption holds.
- **Discrete**: Indicator variables for each post-period (DID-style). Use when effects may not follow a linear trend.
""")

    with st.expander("Analysis Modes"):
        st.markdown("""
- **Calculate MDE given sample size**: You provide cluster counts → app computes the minimum detectable effect size at each power level.
- **Calculate required clusters given target MDE**: You provide a target MDE → app computes how many clusters you need.

**MDE** is expressed in standard deviation units (Cohen's d):

| MDE | Interpretation |
|-----|---------------|
| 0.10 | Very small — requires very large samples |
| 0.20 | Small — typical target for well-powered studies |
| 0.30 | Small-to-medium |
| 0.50 | Medium |
""")

    with st.expander("Timing Groups & Staggered Treatment"):
        st.markdown("""
A **timing group** is a set of clusters that begin treatment at the same time.

- **1 group**: All treated clusters start simultaneously (simple design).
- **2+ groups**: Staggered rollout — different clusters start treatment in different periods.

For each group, set the **start period** (must be ≥ 2, in increasing order) and either:
- **MDE mode**: Number of treatment/comparison clusters in this group.
- **Required clusters mode**: Share (proportion) of total clusters allocated to this group (shares must sum to 1.0).

**Post-period options:**
- *Average across all post-periods*: Usually provides the most power.
- *Single post-period*: Effect at a specific time point or exposure duration.
""")

    with st.expander("Error Structure"):
        st.markdown("""
**Autocorrelation** controls how outcomes within a cluster correlate over time:
- **AR(1)**: Correlation decays geometrically — rho at lag 1, rho² at lag 2, etc.
- **Constant**: Equal correlation (rho) at all lags.
- **None**: No temporal correlation.

**Key parameters:**
- **Cluster autocorrelation (rho)**: Correlation of cluster means over time. Higher rho → *more* variance in DID → *lower* power.
- **Individual autocorrelation (phi)**: Within-person correlation (longitudinal panels only). Higher phi → *less* residual variance → *better* power.
- **ICC**: Proportion of total variance between clusters. Higher ICC → cluster-level variance dominates → need more clusters.
""")

    with st.expander("Precision & Testing"):
        st.markdown("""
- **R²yx** (covariates): Proportion of outcome variance explained by covariates. Higher → less residual variance → better power. Free power if you have good covariates.
- **R²tx** (treatment-covariate correlation): If covariates predict treatment, this *reduces* effective treatment variation → *hurts* power. Usually small in randomized designs.
- **Design effect**: Multiplicative adjustment for survey weights. Default 1.0.
- **Alpha**: Significance level (default 0.05).
- **Two-tailed**: Use for most applications. One-tailed is more powerful but only appropriate with a strong directional hypothesis.
""")

    with st.expander("Scenario Comparison"):
        st.markdown("""
The **Scenario Comparison** tab generates a grid of MDE (or required clusters) across different ICC and R²yx values, holding all other parameters fixed. This is useful for:

1. **Grant proposals**: Show reviewers that your study is well-powered under a range of plausible assumptions.
2. **Sensitivity analysis**: See which parameters most affect your power.
3. **Study planning**: Identify where your design becomes underpowered.

Enter comma-separated ICC and R²yx values, choose a fixed power level, and click **Generate Scenario Table**. The output is designed for easy copy-paste into reports.
""")

    with st.expander("Tips for Practitioners"):
        st.markdown("""
1. **Start with the scenario comparison.** See how sensitive your MDE is to ICC and R²yx before committing to a sample size.
2. **ICC matters more than you think.** Going from ICC = 0.05 to 0.15 can double required sample size.
3. **Covariates are free power.** Even modest R²yx = 0.20 meaningfully improves power without adding clusters.
4. **More time periods help CITS/ITS more than DID.** For CITS/ITS, more periods directly improve slope estimation.
5. **Longitudinal > Cross-sectional for power** (all else equal), but longitudinal panels are harder to collect.
6. **Averaging post-periods** is usually more powerful than looking at a single time point.
""")

    with st.expander("Statistical Formulas"):
        st.markdown("""
This section shows the mathematical formulas underlying each power calculation.
Notation follows Schochet (2022). Subscript $k$ indexes timing groups ($k = 1, \\ldots, K$).

---

### Minimum Detectable Effect (all designs)

$$\\text{MDE} = \\left(t_{1-\\alpha/2,\\,\\nu} + t_{\\beta,\\,\\nu}\\right) \\sqrt{\\text{Var}}$$

- $\\alpha$ = significance level, $\\nu$ = degrees of freedom, $\\beta = 1 - \\text{power}$
- For a one-tailed test, replace $\\alpha/2$ with $\\alpha$
- $t_{p,\\nu}$ is the $p$-th quantile of the $t$-distribution with $\\nu$ df

---

### Precision Factor (all designs)

$$r^2 = \\frac{1 - R^2_{yx}}{1 - R^2_{tx}}$$

- $R^2_{yx}$: proportion of outcome variance explained by covariates → reduces variance
- $R^2_{tx}$: correlation between treatment assignment and covariates → increases variance
- Default: $r^2 = 1$ (no covariates)

---

### DID Variance — Cross-Sectional, Average Post-Period

Let $A_k$ = number of post-periods and $B_k$ = number of pre-periods for group $k$.

**Autocorrelation summary terms** (cluster level, computed from $\\rho$):

$$\\rho_{\\text{post}} = \\frac{2}{A(A-1)} \\sum_{i < j,\\, i,j \\in \\text{post}} \\rho^{|t_j - t_i|}$$

$$\\rho_{\\text{pre}} = \\frac{2}{B(B-1)} \\sum_{i < j,\\, i,j \\in \\text{pre}} \\rho^{|t_j - t_i|}$$

$$\\rho_{\\text{pp}} = \\frac{1}{AB} \\sum_{i \\in \\text{post}} \\sum_{j \\in \\text{pre}} \\rho^{|t_i - t_j|}$$

For a **constant** autocorrelation structure, $\\rho^{|\\text{lag}|}$ is replaced by the flat value $\\rho$.

**Cluster-level variance term:**

$$\\theta_k = \\text{ICC} \\left[ \\frac{1 + \\rho_{\\text{post}}(A_k - 1)}{A_k} + \\frac{1 + \\rho_{\\text{pre}}(B_k - 1)}{B_k} - 2\\rho_{\\text{pp}} \\right]$$

**Individual-level variance term:**

$$e_k = \\frac{1 - \\text{ICC}}{N} \\left( \\frac{1}{A_k} + \\frac{1}{B_k} \\right)$$

**Variance contribution of group $k$:**

$$V_k = A_k^2 \\left(\\frac{1}{M_T} + \\frac{1}{M_C}\\right) r^2 \\left(\\theta_k + e_k\\right)$$

**Total DID variance** (aggregated across $K$ groups):

$$\\text{Var}_{\\text{DID}} = \\text{deff} \\cdot \\frac{\\displaystyle\\sum_k V_k}{\\left(\\displaystyle\\sum_k A_k\\right)^2}$$

---

### Longitudinal Panel Extension

For longitudinal panels, the individual-level term $e_k$ is replaced by a term with the same
structure as $\\theta_k$ but using the **within-person autocorrelation** $\\phi$ in place of $\\rho$,
scaled by $(1-\\text{ICC})/N$ instead of ICC.

---

### CITS Addition — Fully Interacted Specification

CITS adds a trend-correction term on top of the DID formula. Let
$\\Delta_k = \\bar{t}_{\\text{post},k} - \\bar{t}_{\\text{pre},k}$ (difference of period means) and
$S_k = \\sum_{t \\in \\text{pre}} (t - \\bar{t}_{\\text{pre},k})^2$ (sum of squared pre-period time deviations).

The additional cluster-level term is:

$$\\theta_{\\text{CITS},k} = \\text{ICC}\\left[ \\Delta_k^2\\left(\\frac{1}{S_k} + \\frac{(B_k-1)B_k\\,\\tilde{\\rho}_{\\text{pre}}}{S_k^2}\\right) + \\frac{2\\Delta_k B_k\\,\\tilde{\\rho}_{\\text{cross}}}{S_k} - \\frac{2\\Delta_k B_k\\,\\tilde{\\rho}_{\\text{pp}}}{S_k} \\right]$$

where $\\tilde{\\rho}_{\\text{pre}},\\tilde{\\rho}_{\\text{cross}},\\tilde{\\rho}_{\\text{pp}}$ are time-weighted autocorrelation sums (see Schochet 2022, eq. 3–9).

The individual-level addition is:

$$e_{\\text{CITS},k} = \\frac{1-\\text{ICC}}{N} \\cdot \\frac{\\Delta_k^2}{S_k}$$

**ITS** uses the same additions but replaces $(1/M_T + 1/M_C)$ with $1/M_T$ throughout (no comparison group).

---

### Degrees of Freedom

| Design | Specification | Post-period | $\\nu$ |
|--------|--------------|-------------|--------|
| DID | — | Average | $MP - M - KP - \\sum A_k$ |
| DID | — | Specific | $MP - M - n_q P - n_q$ |
| CITS | Fully interacted | Average | $MP - 8K$ |
| CITS | Common slopes | Average | $MP - 6K$ |
| CITS | Discrete | Average | $MP - 4K - \\sum A_k$ |
| ITS | Fully interacted | Average | $MP - 4K$ |
| ITS | Common slopes | Average | $MP - 3K$ |
| ITS | Discrete | Average | $MP - 2K - \\sum A_k$ |

$M$ = total clusters, $P$ = time periods, $K$ = timing groups,
$n_q$ = number of groups for which the specific post-period falls within their post-period range.

---

### Required Clusters (Mode 2)

The app solves for the total $M$ such that $\\text{MDE}(M) = \\text{target MDE}$
using a **secant method** (max 25 iterations, tolerance $10^{-6}$).
Because the degrees of freedom also depend on $M$, the problem is solved iteratively
rather than in closed form.
""")

    st.markdown("---")
    st.markdown(
        "*Based on Schochet, P.Z. (2022). Statistical Power for Estimating Treatment Effects "
        "Using Difference-in-Differences and Comparative Interrupted Time Series Estimators "
        "with Variation in Treatment Timing. "
        "Journal of Educational and Behavioral Statistics, 48(6), 713–751.*"
    )
