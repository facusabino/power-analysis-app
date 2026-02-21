"""
Input validation for the Power Panel app.
Translates the R code's error checking (lines 955-1416).
"""

from dataclasses import dataclass, field
from typing import List
from config import PowerInputs, DesignType, CITSSpec, AnalysisMode


@dataclass
class ValidationResult:
    is_valid: bool = True
    errors: List[str] = field(default_factory=list)

    def add_error(self, msg):
        self.errors.append(msg)
        self.is_valid = False


def validate_inputs(inputs: PowerInputs) -> ValidationResult:
    result = ValidationResult()
    P = inputs.n_time_periods
    K = inputs.n_timing_groups
    sk = inputs.start_times
    did_cits = int(inputs.design_type)
    type_cits = int(inputs.cits_spec)

    # --- Time periods ---
    if P is None or P < 2:
        result.add_error("Number of time periods must be at least 2.")
    elif did_cits > 1 and type_cits in (1, 2) and P < 6:
        result.add_error("Number of time periods must be at least 6 for CITS/ITS with slopes.")
    elif did_cits > 1 and type_cits == 3 and P < 4:
        result.add_error("Number of time periods must be at least 4 for CITS/ITS discrete designs.")

    # --- Timing groups ---
    if K is None or K < 1:
        result.add_error("Number of timing groups must be at least 1.")
    elif len(sk) != K:
        result.add_error(f"Number of start times ({len(sk)}) does not match timing groups ({K}).")
    else:
        prev = 0
        for i, s in enumerate(sk):
            if s is None or s < 2 or (P and s > P):
                result.add_error(f"Invalid start time for timing group {i+1}.")
            elif s <= prev:
                result.add_error(f"Start time for group {i+1} must be after group {i}.")
            else:
                if did_cits > 1 and type_cits < 3 and s < 4:
                    result.add_error("Start times must be at least 4 for CITS/ITS designs.")
                if did_cits > 1 and type_cits < 3 and P and (P - s + 1) < 3:
                    result.add_error(f"Timing group {i+1} must have at least 3 post-periods for CITS/ITS with trendlines.")
            prev = s

    # --- Sample sizes (mode 1) ---
    if inputs.analysis_mode == AnalysisMode.CALC_MDE:
        if inputs.mt is None or inputs.mt < 2:
            result.add_error("Total treatment clusters must be at least 2.")
        if did_cits < 3 and (inputs.mc is None or inputs.mc < 2):
            result.add_error("Total comparison clusters must be at least 2.")

        if len(inputs.mtk) != K:
            result.add_error(f"Treatment clusters per group ({len(inputs.mtk)}) does not match timing groups ({K}).")
        elif inputs.mt and inputs.mt >= 2:
            s = 0
            for i, m in enumerate(inputs.mtk):
                if m is None or m < 2:
                    result.add_error(f"Treatment clusters in group {i+1} must be at least 2.")
                else:
                    s += m
            if s != inputs.mt:
                result.add_error(f"Treatment clusters per group sum ({s}) does not equal total ({inputs.mt}).")

        if did_cits < 3:
            if len(inputs.mck) != K:
                result.add_error(f"Comparison clusters per group ({len(inputs.mck)}) does not match timing groups ({K}).")
            elif inputs.mc and inputs.mc >= 2:
                s = 0
                for i, m in enumerate(inputs.mck):
                    if m is None or m < 2:
                        result.add_error(f"Comparison clusters in group {i+1} must be at least 2.")
                    else:
                        s += m
                if s != inputs.mc:
                    result.add_error(f"Comparison clusters per group sum ({s}) does not equal total ({inputs.mc}).")

    # --- Proportions (mode 2) ---
    if inputs.analysis_mode == AnalysisMode.CALC_CLUSTERS:
        if inputs.target_mde is None or inputs.target_mde <= 0:
            result.add_error("Target MDE must be positive.")
        if did_cits < 3 and (inputs.rt is None or inputs.rt <= 0 or inputs.rt >= 1):
            result.add_error("Treatment proportion must be between 0 and 1 (exclusive).")

        if len(inputs.rtk) != K:
            result.add_error(f"Treatment shares per group ({len(inputs.rtk)}) does not match timing groups ({K}).")
        else:
            s = sum(inputs.rtk)
            if abs(s - 1.0) > 0.01:
                result.add_error(f"Treatment shares must sum to 1 (got {s:.2f}).")

        if did_cits < 3:
            if len(inputs.rck) != K:
                result.add_error(f"Comparison shares per group ({len(inputs.rck)}) does not match timing groups ({K}).")
            else:
                s = sum(inputs.rck)
                if abs(s - 1.0) > 0.01:
                    result.add_error(f"Comparison shares must sum to 1 (got {s:.2f}).")

    # --- Individuals per cluster ---
    if inputs.n_per_cluster is None or inputs.n_per_cluster < 1:
        result.add_error("Individuals per cluster must be at least 1.")

    # --- Scalar parameters ---
    if inputs.icc is None or inputs.icc < 0 or inputs.icc > 1:
        result.add_error("ICC must be between 0 and 1.")
    if inputs.rho is None or inputs.rho <= -1 or inputs.rho >= 1:
        result.add_error("Cluster autocorrelation (rho) must be between -1 and 1 (exclusive).")
    if inputs.phi is None or inputs.phi <= -1 or inputs.phi >= 1:
        result.add_error("Individual autocorrelation (phi) must be between -1 and 1 (exclusive).")
    if inputs.alpha is None or inputs.alpha <= 0 or inputs.alpha >= 1:
        result.add_error("Significance level (alpha) must be between 0 and 1 (exclusive).")
    if inputs.r2yx is None or inputs.r2yx < 0 or inputs.r2yx >= 1:
        result.add_error("R2yx must be between 0 and 1 (exclusive of 1).")
    if inputs.r2tx is None or inputs.r2tx < 0 or inputs.r2tx >= 1:
        result.add_error("R2tx must be between 0 and 1 (exclusive of 1).")
    if inputs.deff_wgt is None or inputs.deff_wgt < 1:
        result.add_error("Design effect must be at least 1.")

    # --- Time intervals (if not evenly spaced) ---
    if not inputs.evenly_spaced:
        ti = inputs.time_intervals
        if ti is None or len(ti) != P:
            result.add_error(f"Time intervals ({len(ti) if ti else 0}) must match number of periods ({P}).")
        elif ti:
            prev_t = 0
            for i, t in enumerate(ti):
                if t is None or t < 1:
                    result.add_error(f"Time interval for period {i+1} must be a positive integer.")
                elif t <= prev_t:
                    result.add_error(f"Time interval for period {i+1} must be after the previous one.")
                prev_t = t

    # --- Post-period specific checks ---
    if inputs.post_period_type.value == 2:
        from config import SpecificPostUnit
        if inputs.specific_post_unit == SpecificPostUnit.TIME_PERIOD:
            q = inputs.q_time
            if q is None or q < 1 or (P and q > P):
                result.add_error("Invalid post-period time point.")
            elif len(sk) > 0 and q < sk[0]:
                result.add_error("Post-period time point must be a post-period for at least one timing group.")
        else:
            l = inputs.l_time
            if l is None or l < 1:
                result.add_error("Invalid exposure time point.")

    return result
