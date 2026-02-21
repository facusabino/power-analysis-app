"""
Independent verification of Power Panel computations.

Two approaches:
1. Benchmark against Schochet (2022) Table 3 values
2. Independent implementation from paper formulas (DID cross-sectional, no AR1)
"""

import numpy as np
from scipy import stats

# =====================================================================
# PART 1: Benchmark against Schochet (2022) Table 3
# =====================================================================
# Table 3 parameters:
#   P=8, K=2, S1=4, S2=6, cross-sectional, AR1, rho=0.4
#   ICC=0.05, N=100, 50-50 split, no covariates, MDE=0.20
#   alpha=0.05, power=0.80, two-tailed
#   Pooled (average across all post-periods)
#
# Paper reports required total clusters (M):
#   DID:                    M = 37
#   CITS fully-interacted:  M = 297
#   CITS common-slopes:     M = 89

print("=" * 70)
print("PART 1: Benchmark against Schochet (2022) Table 3")
print("=" * 70)
print()
print("Paper parameters: P=8, K=2, S1=4, S2=6, cross-sectional, AR1")
print("  rho=0.4, ICC=0.05, N=100, 50-50 split, no covariates")
print("  MDE target=0.20, alpha=0.05, power=0.80, two-tailed")
print()

from config import (
    PowerInputs, DesignType, CITSSpec, AnalysisMode,
    PanelType, AutocorrStructure, PostPeriodType,
)
from power_engine import run_power_analysis

# Common benchmark parameters
benchmark_params = dict(
    n_time_periods=8,
    n_timing_groups=2,
    start_times=[4, 6],
    panel_type=PanelType.CROSS_SECTIONAL,
    autocorr_structure=AutocorrStructure.AR1,
    rho=0.40,
    phi=0.40,
    icc=0.05,
    n_per_cluster=100,
    r2yx=0.0,
    r2tx=0.0,
    deff_wgt=1.0,
    alpha=0.05,
    two_tailed=True,
    evenly_spaced=True,
    post_period_type=PostPeriodType.AVERAGE,
    # Required clusters mode
    analysis_mode=AnalysisMode.CALC_CLUSTERS,
    target_mde=0.20,
    rt=0.50,
    rtk=[0.50, 0.50],
    rck=[0.50, 0.50],
    power_min=0.80,
    power_max=0.80,
    power_step=0.05,
)

# Test DID
inputs_did = PowerInputs(design_type=DesignType.DID, cits_spec=CITSSpec.FULLY_INTERACTED, **benchmark_params)
result_did = run_power_analysis(inputs_did)
m_did = result_did['results'][0]['required_clusters'] if result_did['results'] else None

# Test CITS fully-interacted
inputs_cits_fi = PowerInputs(design_type=DesignType.CITS, cits_spec=CITSSpec.FULLY_INTERACTED, **benchmark_params)
result_cits_fi = run_power_analysis(inputs_cits_fi)
m_cits_fi = result_cits_fi['results'][0]['required_clusters'] if result_cits_fi['results'] else None

# Test CITS common-slopes
inputs_cits_cs = PowerInputs(design_type=DesignType.CITS, cits_spec=CITSSpec.COMMON_SLOPES, **benchmark_params)
result_cits_cs = run_power_analysis(inputs_cits_cs)
m_cits_cs = result_cits_cs['results'][0]['required_clusters'] if result_cits_cs['results'] else None

print(f"{'Design':<30} {'Paper M':>10} {'Our M':>10} {'Match?':>10}")
print("-" * 62)
print(f"{'DID pooled':<30} {'37':>10} {str(m_did):>10} {'OK' if m_did == 37 else 'DIFF':>10}")
print(f"{'CITS fully-interacted':<30} {'297':>10} {str(m_cits_fi):>10} {'OK' if m_cits_fi == 297 else 'DIFF':>10}")
print(f"{'CITS common-slopes':<30} {'89':>10} {str(m_cits_cs):>10} {'OK' if m_cits_cs == 89 else 'DIFF':>10}")
print()

# Additional benchmark: P=12, S1=4, S2=8 from the paper
# DID: M=32, CITS FI: M=641, CITS CS: M=68
print("Additional benchmark: P=12, S1=4, S2=8")
benchmark_params2 = {**benchmark_params, 'n_time_periods': 12, 'start_times': [4, 8]}
inputs_did2 = PowerInputs(design_type=DesignType.DID, cits_spec=CITSSpec.FULLY_INTERACTED, **benchmark_params2)
result_did2 = run_power_analysis(inputs_did2)
m_did2 = result_did2['results'][0]['required_clusters'] if result_did2['results'] else None

inputs_cits_fi2 = PowerInputs(design_type=DesignType.CITS, cits_spec=CITSSpec.FULLY_INTERACTED, **benchmark_params2)
result_cits_fi2 = run_power_analysis(inputs_cits_fi2)
m_cits_fi2 = result_cits_fi2['results'][0]['required_clusters'] if result_cits_fi2['results'] else None

inputs_cits_cs2 = PowerInputs(design_type=DesignType.CITS, cits_spec=CITSSpec.COMMON_SLOPES, **benchmark_params2)
result_cits_cs2 = run_power_analysis(inputs_cits_cs2)
m_cits_cs2 = result_cits_cs2['results'][0]['required_clusters'] if result_cits_cs2['results'] else None

print(f"{'Design':<30} {'Paper M':>10} {'Our M':>10} {'Match?':>10}")
print("-" * 62)
print(f"{'DID pooled':<30} {'32':>10} {str(m_did2):>10} {'OK' if m_did2 == 32 else 'DIFF':>10}")
print(f"{'CITS fully-interacted':<30} {'641':>10} {str(m_cits_fi2):>10} {'OK' if m_cits_fi2 == 641 else 'DIFF':>10}")
print(f"{'CITS common-slopes':<30} {'68':>10} {str(m_cits_cs2):>10} {'OK' if m_cits_cs2 == 68 else 'DIFF':>10}")


# =====================================================================
# PART 2: Independent implementation from paper formulas
# =====================================================================
# We implement the DID cross-sectional variance formula from scratch
# using Equation (13) from the paper for the NO autocorrelation case,
# which reduces to a clean closed form:
#
# For K=1, no AR1 (rho=0), average post-period:
#   Var(beta_DID) = (1/MT + 1/MC) * [ICC*(1/A + 1/B) + (1-ICC)/N * (1/A + 1/B)]
#                 = (1/MT + 1/MC) * (1/A + 1/B) * [ICC + (1-ICC)/N]
#
# MDE = Factor * sqrt(Var)  where Factor = t_inv(0.975, df) + t_inv(0.80, df)

print()
print()
print("=" * 70)
print("PART 2: Independent implementation (DID, no autocorrelation)")
print("=" * 70)
print()
print("Simplified DID formula for K=1, no AR1, cross-sectional, avg post:")
print("  Var = (1/MT + 1/MC) * (1/A + 1/B) * [ICC + (1-ICC)/N]")
print()

def independent_mde_simple(MT, MC, P, S, ICC, N, alpha=0.05, power=0.80):
    """
    Independent MDE calculation using simplified formula.
    K=1 timing group, no autocorrelation, cross-sectional, average post,
    no covariates, even spacing.
    """
    B = S - 1       # pre-periods
    A = P - S + 1   # post-periods
    M = MT + MC
    df = M * P - M - P - A   # DID df for K=1

    var = (1/MT + 1/MC) * (1/A + 1/B) * (ICC + (1 - ICC) / N)

    t_alpha = stats.t.ppf(1 - alpha/2, df)
    t_power = stats.t.ppf(power, df)
    factor = t_alpha + t_power

    mde = factor * np.sqrt(var)
    return mde, var, df

# Test cases with K=1 (no staggered timing)
test_cases = [
    {"MT": 20, "MC": 20, "P": 10, "S": 5, "ICC": 0.05, "N": 100},
    {"MT": 30, "MC": 30, "P": 8, "S": 4, "ICC": 0.10, "N": 200},
    {"MT": 50, "MC": 50, "P": 12, "S": 6, "ICC": 0.02, "N": 500},
    {"MT": 10, "MC": 10, "P": 6, "S": 3, "ICC": 0.15, "N": 50},
]

print(f"{'Case':<8} {'Indep MDE':>12} {'Engine MDE':>12} {'Diff':>10} {'Match?':>8}")
print("-" * 52)

for i, tc in enumerate(test_cases):
    # Independent calculation
    mde_indep, var_indep, df_indep = independent_mde_simple(**tc)

    # Our engine
    inp = PowerInputs(
        design_type=DesignType.DID,
        analysis_mode=AnalysisMode.CALC_MDE,
        panel_type=PanelType.CROSS_SECTIONAL,
        autocorr_structure=AutocorrStructure.NONE,
        n_time_periods=tc["P"],
        n_timing_groups=1,
        start_times=[tc["S"]],
        mt=tc["MT"],
        mtk=[tc["MT"]],
        mc=tc["MC"],
        mck=[tc["MC"]],
        icc=tc["ICC"],
        n_per_cluster=tc["N"],
        rho=0.0, phi=0.0,
        r2yx=0.0, r2tx=0.0,
        deff_wgt=1.0,
        alpha=0.05,
        two_tailed=True,
        evenly_spaced=True,
        post_period_type=PostPeriodType.AVERAGE,
        power_min=0.80,
        power_max=0.80,
        power_step=0.05,
    )
    result = run_power_analysis(inp)
    mde_engine = result['results'][0]['mde_value'] if result['results'] else None

    diff = abs(mde_indep - mde_engine) if mde_engine else float('inf')
    match = "OK" if diff < 0.001 else "DIFF"

    print(f"  {i+1:<6} {mde_indep:>12.4f} {mde_engine:>12.4f} {diff:>10.6f} {match:>8}")


# =====================================================================
# PART 3: Independent implementation with AR1
# =====================================================================
# For K=1, AR1, cross-sectional, average post-period (Equation 13):
#
# Var = (1/MT + 1/MC) * [ICC * theta_term + (1-ICC)/N * e_term]
# where:
#   theta_term = (1 + rho_post*(A-1))/A + (1 + rho_pre*(B-1))/B - 2*rho_pp
#   e_term = 1/A + 1/B
#   rho_pre = avg correlation in pre-period
#   rho_post = avg correlation in post-period
#   rho_pp = avg cross pre-post correlation

print()
print()
print("=" * 70)
print("PART 3: Independent implementation (DID with AR1)")
print("=" * 70)
print()

def calc_rho_independent(rho, indices):
    """Average pairwise AR1 correlation among given 1-based time indices."""
    n = len(indices)
    if n <= 1:
        return 0.0
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            total += rho ** abs(indices[j] - indices[i])
            count += 1
    return total / count if count > 0 else 0.0

def calc_rho_cross(rho, pre_indices, post_indices):
    """Average AR1 correlation across pre and post period indices."""
    total = 0.0
    for i in pre_indices:
        for j in post_indices:
            total += rho ** abs(j - i)
    return total / (len(pre_indices) * len(post_indices))

def independent_mde_ar1(MT, MC, P, S, ICC, N, rho, alpha=0.05, power=0.80):
    """
    Independent MDE: K=1, AR1, cross-sectional, average post, no covariates.
    """
    B = S - 1
    A = P - S + 1
    M = MT + MC
    df = M * P - M - P - A

    pre_idx = list(range(1, S))        # 1..S-1
    post_idx = list(range(S, P + 1))   # S..P

    rho_pre = calc_rho_independent(rho, pre_idx)
    rho_post = calc_rho_independent(rho, post_idx)
    rho_pp = calc_rho_cross(rho, pre_idx, post_idx)

    theta_t1 = (1 + rho_post * (A - 1)) / A
    theta_t2 = (1 + rho_pre * (B - 1)) / B
    theta_t3 = rho_pp
    theta_term = ICC * (theta_t1 + theta_t2 - 2 * theta_t3)

    e_term = ((1 - ICC) / N) * (1.0/A + 1.0/B)

    var = (1/MT + 1/MC) * (theta_term + e_term)

    t_alpha = stats.t.ppf(1 - alpha/2, df)
    t_power = stats.t.ppf(power, df)
    factor = t_alpha + t_power

    mde = factor * np.sqrt(var)
    return mde, var, df

# Test cases with AR1
test_cases_ar1 = [
    {"MT": 20, "MC": 20, "P": 10, "S": 5, "ICC": 0.05, "N": 100, "rho": 0.40},
    {"MT": 30, "MC": 30, "P": 8, "S": 4, "ICC": 0.10, "N": 200, "rho": 0.30},
    {"MT": 50, "MC": 50, "P": 12, "S": 6, "ICC": 0.02, "N": 500, "rho": 0.50},
    {"MT": 15, "MC": 15, "P": 6, "S": 3, "ICC": 0.08, "N": 100, "rho": 0.20},
]

print(f"{'Case':<8} {'Indep MDE':>12} {'Engine MDE':>12} {'Diff':>10} {'Match?':>8}")
print("-" * 52)

for i, tc in enumerate(test_cases_ar1):
    mde_indep, _, _ = independent_mde_ar1(**tc)

    inp = PowerInputs(
        design_type=DesignType.DID,
        analysis_mode=AnalysisMode.CALC_MDE,
        panel_type=PanelType.CROSS_SECTIONAL,
        autocorr_structure=AutocorrStructure.AR1,
        n_time_periods=tc["P"],
        n_timing_groups=1,
        start_times=[tc["S"]],
        mt=tc["MT"],
        mtk=[tc["MT"]],
        mc=tc["MC"],
        mck=[tc["MC"]],
        icc=tc["ICC"],
        n_per_cluster=tc["N"],
        rho=tc["rho"], phi=0.0,
        r2yx=0.0, r2tx=0.0,
        deff_wgt=1.0,
        alpha=0.05,
        two_tailed=True,
        evenly_spaced=True,
        post_period_type=PostPeriodType.AVERAGE,
        power_min=0.80,
        power_max=0.80,
        power_step=0.05,
    )
    result = run_power_analysis(inp)
    mde_engine = result['results'][0]['mde_value'] if result['results'] else None

    diff = abs(mde_indep - mde_engine) if mde_engine else float('inf')
    match = "OK" if diff < 0.001 else "DIFF"

    print(f"  {i+1:<6} {mde_indep:>12.4f} {mde_engine:>12.4f} {diff:>10.6f} {match:>8}")


print()
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("Part 1: Benchmarks against Schochet (2022) Table 3")
print("  - Compares required cluster counts for DID, CITS FI, CITS CS")
print("  - Using the paper's exact parameter configurations")
print()
print("Part 2: Independent formula (no AR1) vs engine")
print("  - Implements Var formula from scratch (Equation 13, simplified)")
print("  - Tests across 4 different parameter sets")
print()
print("Part 3: Independent formula (with AR1) vs engine")
print("  - Implements full AR1 autocorrelation from scratch")
print("  - Tests across 4 different parameter sets")
