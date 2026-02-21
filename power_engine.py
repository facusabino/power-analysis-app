"""
Core power analysis computation module for panel data designs.
Translates the R Shiny app logic from Schochet (JEBS, 2022).
No Streamlit dependency — pure computation.
"""

import numpy as np
from scipy import stats
from config import (
    PowerInputs, DesignType, CITSSpec, AnalysisMode,
    PanelType, AutocorrStructure, PostPeriodType, SpecificPostUnit,
)


def calc_rho(corrp, prepost, pres, pref, posts, postf, tbarb1, tbarb2, ar, tpf):
    """
    Calculate autocorrelation and time terms for variance/covariance components.

    Parameters
    ----------
    corrp : float – correlation parameter (rho or phi)
    prepost : int – case selector (1–12, 90, 100, 110)
    pres, pref : int – pre-period start/end (1-based)
    posts, postf : int – post-period start/end (1-based)
    tbarb1, tbarb2 : float – time-bar means
    ar : int – 1 for AR1, 0 for constant correlation
    tpf : np.ndarray – time period values, 1-based (index 0 is unused padding)
    """

    if prepost == 1:
        # Pre-period within-group autocorrelation
        rho_ret = 0.0
        btemp = pref - pres + 1
        if btemp > 1:
            indices = np.arange(pres, pref + 1)
            pairs = np.array([(i, j) for i in indices for j in indices if j > i])
            if ar == 1:
                sv = corrp ** (tpf[pairs[:, 1]] - tpf[pairs[:, 0]])
            else:
                sv = np.full(len(pairs), corrp)
            rho_ret = 2.0 * np.sum(sv) / (btemp * (btemp - 1))
        return rho_ret

    elif prepost == 2:
        # Post-period within-group autocorrelation
        rho_ret = 0.0
        atemp = postf - posts + 1
        if atemp > 1:
            indices = np.arange(posts, postf + 1)
            pairs = np.array([(i, j) for i in indices for j in indices if j > i])
            if ar == 1:
                sv = corrp ** (tpf[pairs[:, 1]] - tpf[pairs[:, 0]])
            else:
                sv = np.full(len(pairs), corrp)
            rho_ret = 2.0 * np.sum(sv) / (atemp * (atemp - 1))
        return rho_ret

    elif prepost == 3:
        # Cross pre-post autocorrelation
        btemp = pref - pres + 1
        atemp = postf - posts + 1
        pre_idx = np.arange(pres, pref + 1)
        post_idx = np.arange(posts, postf + 1)
        grid = np.array([(i, j) for i in pre_idx for j in post_idx])
        if ar == 1:
            sv = corrp ** (tpf[grid[:, 1]] - tpf[grid[:, 0]])
        else:
            sv = np.full(len(grid), corrp)
        return np.sum(sv) / (atemp * btemp)

    elif prepost == 4:
        # Post-period time difference sum
        rho_ret = 0.0
        atemp = postf - posts + 1
        if atemp > 1:
            indices = np.arange(posts, postf + 1)
            sv = tpf[indices] - tbarb1
            rho_ret = np.sum(sv)
        return rho_ret

    elif prepost == 5:
        # Post-period covariance product sum
        rho_ret = 0.0
        atemp = postf - posts + 1
        if atemp >= 1:
            indices = np.arange(posts, postf + 1)
            grid = np.array([(i, j) for i in indices for j in indices])
            sv = (tpf[grid[:, 0]] - tbarb1) * (tpf[grid[:, 1]] - tbarb1)
            rho_ret = np.sum(sv)
        return rho_ret

    elif prepost == 6:
        # Pre-period time-weighted autocorrelation (for CITS)
        rho_ret = 0.0
        btemp = pref - pres + 1
        if btemp > 1:
            indices = np.arange(pres, pref + 1)
            pairs = np.array([(i, j) for i in indices for j in indices if j > i])
            time_term = (tpf[pairs[:, 0]] - tbarb1) * (tpf[pairs[:, 1]] - tbarb1)
            if ar == 1:
                sv = time_term * (corrp ** (tpf[pairs[:, 1]] - tpf[pairs[:, 0]]))
            else:
                sv = time_term * corrp
            rho_ret = 2.0 * np.sum(sv) / (btemp * (btemp - 1))
        return rho_ret

    elif prepost == 7:
        # Pre-period single time-weight autocorrelation
        rho_ret = 0.0
        btemp = pref - pres + 1
        if btemp > 1:
            indices = np.arange(pres, pref + 1)
            grid = np.array([(i, j) for i in indices for j in indices])
            time_term = tpf[grid[:, 0]] - tbarb1
            if ar == 1:
                sv = time_term * (corrp ** np.abs(tpf[grid[:, 1]] - tpf[grid[:, 0]]))
            else:
                sv = time_term * corrp
            rho_ret = np.sum(sv) / (btemp * btemp)
        return rho_ret

    elif prepost == 8:
        # Cross pre-post time-weighted
        btemp = pref - pres + 1
        atemp = postf - posts + 1
        pre_idx = np.arange(pres, pref + 1)
        post_idx = np.arange(posts, postf + 1)
        grid = np.array([(i, j) for i in pre_idx for j in post_idx])
        time_term = tpf[grid[:, 0]] - tbarb1
        if ar == 1:
            sv = time_term * (corrp ** (tpf[grid[:, 1]] - tpf[grid[:, 0]]))
        else:
            sv = time_term * corrp
        return np.sum(sv) / (btemp * atemp)

    elif prepost == 9:
        # Post-period time-weighted autocorrelation
        rho_ret = 0.0
        atemp = postf - posts + 1
        if atemp > 1:
            indices = np.arange(posts, postf + 1)
            pairs = np.array([(i, j) for i in indices for j in indices if j > i])
            time_term = (tpf[pairs[:, 0]] - tbarb1) * (tpf[pairs[:, 1]] - tbarb1)
            if ar == 1:
                sv = time_term * (corrp ** (tpf[pairs[:, 1]] - tpf[pairs[:, 0]]))
            else:
                sv = time_term * corrp
            rho_ret = 2.0 * np.sum(sv) / (atemp * (atemp - 1))
        return rho_ret

    elif prepost == 10:
        # Post-period single time-weight
        rho_ret = 0.0
        atemp = postf - posts + 1
        if atemp > 1:
            indices = np.arange(posts, postf + 1)
            grid = np.array([(i, j) for i in indices for j in indices])
            time_term = tpf[grid[:, 0]] - tbarb1
            if ar == 1:
                sv = time_term * (corrp ** np.abs(tpf[grid[:, 1]] - tpf[grid[:, 0]]))
            else:
                sv = time_term * corrp
            rho_ret = np.sum(sv) / (atemp * atemp)
        return rho_ret

    elif prepost == 11:
        # Cross pre-post with post time-weight
        btemp = pref - pres + 1
        atemp = postf - posts + 1
        pre_idx = np.arange(pres, pref + 1)
        post_idx = np.arange(posts, postf + 1)
        grid = np.array([(i, j) for i in pre_idx for j in post_idx])
        time_term = tpf[grid[:, 1]] - tbarb1
        if ar == 1:
            sv = time_term * (corrp ** (tpf[grid[:, 1]] - tpf[grid[:, 0]]))
        else:
            sv = time_term * corrp
        return np.sum(sv) / (btemp * atemp)

    elif prepost == 12:
        # Cross pre-post with both time-weights
        btemp = pref - pres + 1
        atemp = postf - posts + 1
        pre_idx = np.arange(pres, pref + 1)
        post_idx = np.arange(posts, postf + 1)
        grid = np.array([(i, j) for i in pre_idx for j in post_idx])
        time_term = (tpf[grid[:, 0]] - tbarb1) * (tpf[grid[:, 1]] - tbarb2)
        if ar == 1:
            sv = time_term * (corrp ** (tpf[grid[:, 1]] - tpf[grid[:, 0]]))
        else:
            sv = time_term * corrp
        return np.sum(sv) / (btemp * atemp)

    elif prepost == 90:
        # Pre-period product of two different time-bar deviations
        rho_ret = 0.0
        btemp = pref - pres + 1
        if btemp > 1:
            indices = np.arange(pres, pref + 1)
            sv = (tpf[indices] - tbarb1) * (tpf[indices] - tbarb2)
            rho_ret = np.sum(sv)
        return rho_ret

    elif prepost == 100:
        # Pre-period cross time-weight autocorrelation
        rho_ret = 0.0
        btemp = pref - pres + 1
        if btemp > 1:
            indices = np.arange(pres, pref + 1)
            pairs = np.array([(i, j) for i in indices for j in indices if j != i])
            time_term = (tpf[pairs[:, 0]] - tbarb1) * (tpf[pairs[:, 1]] - tbarb2)
            if ar == 1:
                sv = time_term * (corrp ** np.abs(tpf[pairs[:, 1]] - tpf[pairs[:, 0]]))
            else:
                sv = time_term * corrp
            rho_ret = np.sum(sv) / (btemp * (btemp - 1))
        return rho_ret

    elif prepost == 110:
        # Cross pre-post with both different time-bar weights
        btemp = pref - pres + 1
        atemp = postf - posts + 1
        pre_idx = np.arange(pres, pref + 1)
        post_idx = np.arange(posts, postf + 1)
        grid = np.array([(i, j) for i in pre_idx for j in post_idx])
        time_term = (tpf[grid[:, 0]] - tbarb1) * (tpf[grid[:, 1]] - tbarb2)
        if ar == 1:
            sv = time_term * (corrp ** (tpf[grid[:, 1]] - tpf[grid[:, 0]]))
        else:
            sv = time_term * corrp
        return np.sum(sv) / (btemp * atemp)

    return 0.0


def setup_time_periods(inputs):
    """Compute derived time variables for each timing group."""
    P = inputs.n_time_periods
    K = inputs.n_timing_groups
    sk = inputs.start_times

    # Build 1-based time period array (index 0 = padding)
    if inputs.evenly_spaced:
        tpp = np.zeros(P + 1)
        for i in range(1, P + 1):
            tpp[i] = i
    else:
        tpp = np.zeros(P + 1)
        for i in range(P):
            tpp[i + 1] = inputs.time_intervals[i]

    # Per-group quantities
    tbar = np.zeros(K)
    ssqt = np.zeros(K)
    tbara = np.zeros(K)
    ssqta = np.zeros(K)
    pp_diff = np.zeros(K)
    bk = np.zeros(K, dtype=int)
    ak = np.zeros(K, dtype=int)

    for k in range(K):
        stk = sk[k]
        bk[k] = stk - 1
        ak[k] = P - stk + 1

        # Pre-period: indices 1..stk-1
        pre_times = tpp[1:stk]
        tbar[k] = np.mean(pre_times)
        ssqt[k] = np.sum((pre_times - tbar[k]) ** 2)

        # Post-period: indices stk..P
        post_times = tpp[stk:P + 1]
        tbara[k] = np.mean(post_times)
        ssqta[k] = np.sum((post_times - tbara[k]) ** 2)

        pp_diff[k] = tbara[k] - tbar[k]

    # Full-period
    all_times = tpp[1:P + 1]
    tbar_full = np.mean(all_times)
    ssqt_full = np.sum((all_times - tbar_full) ** 2)

    return {
        'tpp': tpp, 'tbar': tbar, 'ssqt': ssqt,
        'tbara': tbara, 'ssqta': ssqta, 'pp_diff': pp_diff,
        'bk': bk, 'ak': ak, 'tbar_full': tbar_full, 'ssqt_full': ssqt_full,
    }


def compute_df(design, cits_spec, m_total, n_time, n_timing_groups, sum_ak, sum_iq, avg_point):
    """Compute degrees of freedom for the t-distribution."""
    M = m_total
    P = n_time
    K = n_timing_groups

    if avg_point == 1:
        n_groups = K
    else:
        n_groups = sum_iq

    if design == DesignType.DID:
        if avg_point == 1:
            df = M * P - M - K * P - sum_ak
        else:
            df = M * P - M - sum_iq * P - sum_iq
    elif design == DesignType.CITS:
        if cits_spec == CITSSpec.FULLY_INTERACTED:
            df = M * P - 8 * n_groups
        elif cits_spec == CITSSpec.COMMON_SLOPES:
            df = M * P - 6 * n_groups
        elif cits_spec == CITSSpec.DISCRETE:
            if avg_point == 1:
                df = M * P - 4 * K - sum_ak
            else:
                df = M * P - 4 * sum_iq - sum_iq
    elif design == DesignType.ITS:
        if cits_spec == CITSSpec.FULLY_INTERACTED:
            df = M * P - 4 * n_groups
        elif cits_spec == CITSSpec.COMMON_SLOPES:
            df = M * P - 3 * n_groups
        elif cits_spec == CITSSpec.DISCRETE:
            if avg_point == 1:
                df = M * P - 2 * K - sum_ak
            else:
                df = M * P - 2 * sum_iq - sum_iq

    return max(df, 1)


def compute_mde(variance, df, alpha, power, two_tailed):
    """MDE = Factor(alpha, power, df) * sqrt(variance)"""
    alpha_adj = 1.0 - (alpha / 2.0) if two_tailed else 1.0 - alpha
    inv_alpha = stats.t.ppf(alpha_adj, df)
    inv_power = stats.t.ppf(power, df)
    factor = inv_alpha + inv_power
    return factor * np.sqrt(variance)


def compute_required_clusters(var_term, target_mde, alpha, power, two_tailed,
                              df_func, max_iter=25, tol=1e-6):
    """
    Use secant method to find total M such that MDE(M) = target_mde.
    df_func(m) returns degrees of freedom as a function of total clusters m.
    """
    if var_term <= 0:
        return {'m_opt': None, 'converged': False}

    alpha_adj = 1.0 - (alpha / 2.0) if two_tailed else 1.0 - alpha
    m0 = 30.0
    m1 = 20.0

    def f(m):
        df = max(df_func(m), 1)
        inv_a = stats.t.ppf(alpha_adj, df)
        inv_p = stats.t.ppf(power, df)
        fac = inv_a + inv_p
        return m - (fac ** 2) * var_term / (target_mde ** 2)

    f0 = f(m0)
    f1 = f(m1)

    for iteration in range(max_iter):
        if abs(f0) <= tol:
            return {'m_opt': round(m0), 'converged': True}
        if abs(f0 - f1) < 1e-15:
            break
        delta = f0 * (m0 - m1) / (f0 - f1)
        m1 = m0
        f1 = f0
        m0 = m0 - delta
        f0 = f(m0)

    if abs(f0) <= tol:
        return {'m_opt': round(m0), 'converged': True}
    return {'m_opt': None, 'converged': False}


def compute_variance(inputs, time_info):
    """
    Main variance calculation loop over timing groups.
    Returns dict with 'did_tot', 'cits_tot', 'its_tot'.
    """
    P = inputs.n_time_periods
    K = inputs.n_timing_groups
    sk = inputs.start_times
    tpp = time_info['tpp']
    bk = time_info['bk']
    ak = time_info['ak']
    tbar = time_info['tbar']
    ssqt_arr = time_info['ssqt']
    tbara = time_info['tbara']
    ssqta = time_info['ssqta']
    pp_diff = time_info['pp_diff']
    tbar_full = time_info['tbar_full']
    ssqt_full = time_info['ssqt_full']

    icc = inputs.icc
    nsamp = inputs.n_per_cluster
    rho = inputs.rho
    phi = inputs.phi
    deff_wgt = inputs.deff_wgt

    did_cits = int(inputs.design_type)
    type_cits = int(inputs.cits_spec)
    cross_sec = int(inputs.panel_type)
    avg_point = int(inputs.post_period_type)

    # Handle autocorrelation structure
    ar1_orig = int(inputs.autocorr_structure)
    if ar1_orig == 3:  # None
        rho = 0.0
        phi = 0.0
    ar1 = 1 if ar1_orig == 1 else 0  # AR1=1 -> ar=1; Constant or None -> ar=0

    # Precision factor
    if 0 <= inputs.r2yx < 1 and 0 <= inputs.r2tx < 1:
        r2 = (1 - inputs.r2yx) / (1 - inputs.r2tx)
    else:
        r2 = 1.0

    # Build cluster arrays
    if inputs.analysis_mode == AnalysisMode.CALC_CLUSTERS:
        mt = inputs.rt
        mtk = [inputs.rt * r for r in inputs.rtk]
        mtk_its = list(inputs.rtk)
        mc = 1 - inputs.rt
        mck = [(1 - inputs.rt) * r for r in inputs.rck]
    else:
        mt = inputs.mt
        mtk = list(inputs.mtk)
        mtk_its = list(inputs.mtk)
        mc = inputs.mc
        mck = list(inputs.mck)

    # Weight vector and sums
    sumak = sum(int(a) for a in ak)
    sumwk = sumak

    # Initialize accumulation terms
    term1 = 0.0
    term1_long = 0.0
    term1_noc = 0.0
    term1_long_noc = 0.0

    term1q = 0.0
    term1q_long = 0.0
    term1q_noc = 0.0
    term1q_long_noc = 0.0

    sumiq = 0
    cits_term1 = 0.0
    its_term1 = 0.0
    cits_term1_long = 0.0
    its_term1_long = 0.0

    q_time = inputs.q_time
    l_time = inputs.l_time

    for k in range(K):
        bkk = int(bk[k])
        skk = sk[k]
        akk = int(ak[k])
        mtkk = mtk[k]
        mckk = mck[k] if did_cits < 3 else 0
        mtkk_its = mtk_its[k]
        tbarkk = tbar[k]
        ssqtkk = ssqt_arr[k]
        tbarakk = tbara[k]
        ssqtakk = ssqta[k]
        pp_diffkk = pp_diff[k]

        iqkk = 1

        # ---- DID for average post-period ----
        if avg_point == 1:
            rho_pre = calc_rho(rho, 1, 1, bkk, 0, 0, 0, 0, ar1, tpp)
            rho_post = calc_rho(rho, 2, 0, 0, skk, P, 0, 0, ar1, tpp)
            rho_pp = calc_rho(rho, 3, 1, bkk, skk, P, 0, 0, ar1, tpp)

            theta_term1 = (1 + rho_post * (akk - 1)) / akk
            theta_term2 = (1 + rho_pre * (bkk - 1)) / bkk
            theta_term3 = rho_pp
            theta_term = icc * (theta_term1 + theta_term2 - 2 * theta_term3)

            e_term = ((1 - icc) / nsamp) * (1.0 / akk + 1.0 / bkk)

            inv_clus = (1.0 / mtkk) + (1.0 / mckk) if did_cits < 3 else (1.0 / mtkk_its)

            term1 += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_term + e_term) if did_cits < 3 else 0
            term1_noc += (akk ** 2) * (1.0 / mtkk_its) * r2 * (theta_term + e_term)

            if did_cits < 3:
                term1 += 0  # already added above
            # Fix: only add to term1 once
            # (The above logic handles DID and ITS separately)

            if cross_sec == 2:
                phi_pre = calc_rho(phi, 1, 1, bkk, 0, 0, 0, 0, ar1, tpp)
                phi_post = calc_rho(phi, 2, 0, 0, skk, P, 0, 0, ar1, tpp)
                phi_pp = calc_rho(phi, 3, 1, bkk, skk, P, 0, 0, ar1, tpp)

                e_term1_l = (1 + phi_post * (akk - 1)) / akk
                e_term2_l = (1 + phi_pre * (bkk - 1)) / bkk
                e_term3_l = phi_pp
                e_term_long = ((1 - icc) / nsamp) * (e_term1_l + e_term2_l - 2 * e_term3_l)

                if did_cits < 3:
                    term1_long += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_term + e_term_long)
                term1_long_noc += (akk ** 2) * (1.0 / mtkk_its) * r2 * (theta_term + e_term_long)

        # ---- DID for specific post-period ----
        if avg_point == 2:
            if inputs.specific_post_unit == SpecificPostUnit.TIME_PERIOD:
                qkk = q_time
            else:
                qkk = l_time + skk - 1

            if skk <= qkk <= P:
                iqkk = 1
                sumiq += 1
            else:
                iqkk = 0

            if iqkk == 1:
                rho_preq = calc_rho(rho, 1, 1, bkk, 0, 0, 0, 0, ar1, tpp)
                rho_ppq = calc_rho(rho, 3, 1, bkk, qkk, qkk, 0, 0, ar1, tpp)

                theta_term1q = 1.0
                theta_term2q = (1 + rho_preq * (bkk - 1)) / bkk
                theta_term3q = rho_ppq
                theta_termq = icc * (theta_term1q + theta_term2q - 2 * theta_term3q)

                e_termq = ((1 - icc) / nsamp) * (1 + 1.0 / bkk)

                if did_cits < 3:
                    term1q += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_termq + e_termq)
                term1q_noc += (1.0 / mtkk_its) * r2 * (theta_termq + e_termq)

                # Also compute the average-based DID variance needed for CITS
                rho_pre = calc_rho(rho, 1, 1, bkk, 0, 0, 0, 0, ar1, tpp)
                rho_post = calc_rho(rho, 2, 0, 0, skk, P, 0, 0, ar1, tpp)
                rho_pp = calc_rho(rho, 3, 1, bkk, skk, P, 0, 0, ar1, tpp)

                theta_term1_a = (1 + rho_post * (akk - 1)) / akk
                theta_term2_a = (1 + rho_pre * (bkk - 1)) / bkk
                theta_term3_a = rho_pp
                theta_term_a = icc * (theta_term1_a + theta_term2_a - 2 * theta_term3_a)
                e_term_a = ((1 - icc) / nsamp) * (1.0 / akk + 1.0 / bkk)

                if did_cits < 3:
                    term1 += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_term_a + e_term_a)
                term1_noc += (1.0 / mtkk_its) * r2 * (theta_term_a + e_term_a)

                if cross_sec == 2:
                    phi_pre = calc_rho(phi, 1, 1, bkk, 0, 0, 0, 0, ar1, tpp)
                    phi_post = calc_rho(phi, 2, 0, 0, skk, P, 0, 0, ar1, tpp)
                    phi_pp_q = calc_rho(phi, 3, 1, bkk, qkk, qkk, 0, 0, ar1, tpp)

                    e_term1q_l = 1.0
                    e_term2q_l = (1 + phi_pre * (bkk - 1)) / bkk
                    e_term3q_l = phi_pp_q
                    e_term_longq = ((1 - icc) / nsamp) * (e_term1q_l + e_term2q_l - 2 * e_term3q_l)

                    if did_cits < 3:
                        term1q_long += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_termq + e_term_longq)
                    term1q_long_noc += (1.0 / mtkk_its) * r2 * (theta_termq + e_term_longq)

                    phi_post_a = calc_rho(phi, 2, 0, 0, skk, P, 0, 0, ar1, tpp)
                    phi_pp_a = calc_rho(phi, 3, 1, bkk, skk, P, 0, 0, ar1, tpp)
                    e_term1_la = (1 + phi_post_a * (akk - 1)) / akk
                    e_term2_la = (1 + phi_pre * (bkk - 1)) / bkk
                    e_term3_la = phi_pp_a
                    e_term_long_a = ((1 - icc) / nsamp) * (e_term1_la + e_term2_la - 2 * e_term3_la)

                    if did_cits < 3:
                        term1_long += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (theta_term_a + e_term_long_a)
                    term1_long_noc += (1.0 / mtkk_its) * r2 * (theta_term_a + e_term_long_a)

        # ---- CITS fully interacted, avg post ----
        if did_cits > 1 and type_cits == 1 and avg_point == 1:
            rho_pre1 = calc_rho(rho, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pre2 = calc_rho(rho, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pp1 = calc_rho(rho, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)

            ct1 = (pp_diffkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * rho_pre1 / (ssqtkk ** 2)))
            ct2 = 2 * pp_diffkk * bkk * rho_pre2 / ssqtkk
            ct3 = 2 * pp_diffkk * bkk * rho_pp1 / ssqtkk
            cits_theta = icc * (ct1 + ct2 - ct3)

            ce4 = (pp_diffkk ** 2) / ssqtkk
            cits_e = ((1 - icc) / nsamp) * ce4

            if did_cits == 2:
                cits_term1 += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
            its_term1 += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)

            if cross_sec == 2:
                phi_pre1 = calc_rho(phi, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pre2 = calc_rho(phi, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pp1 = calc_rho(phi, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)

                ce1 = (pp_diffkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * phi_pre1 / (ssqtkk ** 2)))
                ce2 = 2 * pp_diffkk * bkk * phi_pre2 / ssqtkk
                ce3 = 2 * pp_diffkk * bkk * phi_pp1 / ssqtkk
                cits_e_long = ((1 - icc) / nsamp) * (ce1 + ce2 - ce3)

                if did_cits == 2:
                    cits_term1_long += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                its_term1_long += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)

        # ---- CITS fully interacted, specific post ----
        if did_cits > 1 and type_cits == 1 and avg_point == 2 and iqkk == 1:
            rho_pre1 = calc_rho(rho, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pre2 = calc_rho(rho, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_post1 = calc_rho(rho, 9, 0, 0, skk, P, tbarakk, 0, ar1, tpp)
            rho_post2 = calc_rho(rho, 10, 0, 0, skk, P, tbarakk, 0, ar1, tpp)
            rho_pp2 = calc_rho(rho, 11, 1, bkk, skk, P, tbarakk, 0, ar1, tpp)
            rho_pp3 = calc_rho(rho, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)
            rho_pp4 = calc_rho(rho, 12, 1, bkk, skk, P, tbarkk, tbarakk, ar1, tpp)

            diffa_qkk = tpp[qkk] - tbarakk
            diff_qkk = tpp[qkk] - tbarkk

            ct1 = (diffa_qkk ** 2) * ((1 / ssqtakk) + ((akk - 1) * akk * rho_post1 / (ssqtakk ** 2)))
            ct2 = (diff_qkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * rho_pre1 / (ssqtkk ** 2)))
            ct3 = 2 * diffa_qkk * akk * rho_post2 / ssqtakk
            ct4 = 2 * diff_qkk * bkk * rho_pre2 / ssqtkk
            ct5 = 2 * diffa_qkk * akk * rho_pp2 / ssqtakk
            ct6 = 2 * diff_qkk * bkk * rho_pp3 / ssqtkk
            ct7 = 2 * diffa_qkk * diff_qkk * akk * bkk * rho_pp4 / (ssqtakk * ssqtkk)
            cits_theta = icc * (ct1 + ct2 + ct3 + ct4 - ct5 - ct6 - ct7)

            ce8 = (diffa_qkk ** 2) / ssqtakk
            ce9 = (diff_qkk ** 2) / ssqtkk
            cits_e = ((1 - icc) / nsamp) * (ce8 + ce9)

            if did_cits == 2:
                cits_term1 += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
            its_term1 += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)

            if cross_sec == 2:
                phi_pre1 = calc_rho(phi, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pre2 = calc_rho(phi, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_post1 = calc_rho(phi, 9, 0, 0, skk, P, tbarakk, 0, ar1, tpp)
                phi_post2 = calc_rho(phi, 10, 0, 0, skk, P, tbarakk, 0, ar1, tpp)
                phi_pp2 = calc_rho(phi, 11, 1, bkk, skk, P, tbarakk, 0, ar1, tpp)
                phi_pp3 = calc_rho(phi, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)
                phi_pp4 = calc_rho(phi, 12, 1, bkk, skk, P, tbarkk, tbarakk, ar1, tpp)

                ce1 = (diffa_qkk ** 2) * ((1 / ssqtakk) + ((akk - 1) * akk * phi_post1 / (ssqtakk ** 2)))
                ce2 = (diff_qkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * phi_pre1 / (ssqtkk ** 2)))
                ce3 = 2 * diffa_qkk * akk * phi_post2 / ssqtakk
                ce4 = 2 * diff_qkk * bkk * phi_pre2 / ssqtkk
                ce5 = 2 * diffa_qkk * akk * phi_pp2 / ssqtakk
                ce6 = 2 * diff_qkk * bkk * phi_pp3 / ssqtkk
                ce7 = 2 * diffa_qkk * diff_qkk * akk * bkk * phi_pp4 / (ssqtakk * ssqtkk)
                cits_e_long = ((1 - icc) / nsamp) * (ce1 + ce2 + ce3 + ce4 - ce5 - ce6 - ce7)

                if did_cits == 2:
                    cits_term1_long += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                its_term1_long += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)

        # ---- CITS common slopes ----
        if did_cits > 1 and type_cits == 2:
            postd = np.zeros(P + 1)
            postd[skk:P + 1] = 1.0
            postmean = akk / P

            # rho_full1: correlations of post-indicator deviations
            rho_full1 = 0.0
            if P > 1:
                indices = np.arange(1, P + 1)
                pairs = np.array([(i, j) for i in indices for j in indices if j > i])
                time_term = (postd[pairs[:, 0]] - postmean) * (postd[pairs[:, 1]] - postmean)
                if ar1 == 1:
                    sv = time_term * (rho ** (tpp[pairs[:, 1]] - tpp[pairs[:, 0]]))
                else:
                    sv = time_term * rho
                rho_full1 = 2 * np.sum(sv) / (P * (P - 1))

            # rho_full2: cross time-post correlations
            rho_full2 = 0.0
            if P > 1:
                indices = np.arange(1, P + 1)
                pairs_ne = np.array([(i, j) for i in indices for j in indices if j != i])
                time_term = (tpp[pairs_ne[:, 0]] - tbar_full) * (postd[pairs_ne[:, 1]] - postmean)
                if ar1 == 1:
                    sv = time_term * (rho ** np.abs(tpp[pairs_ne[:, 1]] - tpp[pairs_ne[:, 0]]))
                else:
                    sv = time_term * rho
                rho_full2 = np.sum(sv) / (P * (P - 1))

            rho_full3 = calc_rho(rho, 6, 1, P, 0, 0, tbar_full, 0, ar1, tpp)

            ssqt_term = ssqt_full / (ssqtkk + ssqtakk)
            ct1 = ssqt_term
            ct2 = (1.0 / akk + 1.0 / bkk) * P * (P - 1) * (ssqt_term ** 2) * rho_full1
            ct3 = 2 * P * (P - 1) * (ssqt_full / (ssqtkk + ssqtakk)) * (pp_diffkk / (ssqtkk + ssqtakk)) * rho_full2
            ct4 = akk * bkk * (P - 1) * ((pp_diffkk / (ssqtkk + ssqtakk)) ** 2) * rho_full3
            cits_theta = icc * (1.0 / akk + 1.0 / bkk) * (ct1 + ct2 - ct3 + ct4)

            ce1 = ssqt_term
            cits_e = ((1 - icc) / nsamp) * (1.0 / akk + 1.0 / bkk) * ce1

            if avg_point == 1:
                if did_cits == 2:
                    cits_term1 += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
                its_term1 += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)
            elif avg_point == 2 and iqkk == 1:
                if did_cits == 2:
                    cits_term1 += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
                its_term1 += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)

            # Longitudinal
            if cross_sec == 2:
                phi_full1 = 0.0
                if P > 1:
                    indices = np.arange(1, P + 1)
                    pairs = np.array([(i, j) for i in indices for j in indices if j > i])
                    time_term = (postd[pairs[:, 0]] - postmean) * (postd[pairs[:, 1]] - postmean)
                    if ar1 == 1:
                        sv = time_term * (phi ** (tpp[pairs[:, 1]] - tpp[pairs[:, 0]]))
                    else:
                        sv = time_term * phi
                    phi_full1 = 2 * np.sum(sv) / (P * (P - 1))

                phi_full2 = 0.0
                if P > 1:
                    indices = np.arange(1, P + 1)
                    pairs_ne = np.array([(i, j) for i in indices for j in indices if j != i])
                    time_term = (tpp[pairs_ne[:, 0]] - tbar_full) * (postd[pairs_ne[:, 1]] - postmean)
                    if ar1 == 1:
                        sv = time_term * (phi ** np.abs(tpp[pairs_ne[:, 1]] - tpp[pairs_ne[:, 0]]))
                    else:
                        sv = time_term * phi
                    phi_full2 = np.sum(sv) / (P * (P - 1))

                phi_full3 = calc_rho(phi, 6, 1, P, 0, 0, tbar_full, 0, ar1, tpp)

                ce1_l = ssqt_term
                ce2_l = (1.0 / akk + 1.0 / bkk) * P * (P - 1) * (ssqt_term ** 2) * phi_full1
                ce3_l = 2 * P * (P - 1) * (ssqt_full / (ssqtkk + ssqtakk)) * (pp_diffkk / (ssqtkk + ssqtakk)) * phi_full2
                ce4_l = akk * bkk * (P - 1) * ((pp_diffkk / (ssqtkk + ssqtakk)) ** 2) * phi_full3
                cits_e_long = ((1 - icc) / nsamp) * (1.0 / akk + 1.0 / bkk) * (ce1_l + ce2_l - ce3_l + ce4_l)

                if avg_point == 1:
                    if did_cits == 2:
                        cits_term1_long += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                    its_term1_long += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)
                elif avg_point == 2 and iqkk == 1:
                    if did_cits == 2:
                        cits_term1_long += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                    its_term1_long += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)

        # ---- CITS discrete, avg post ----
        if did_cits > 1 and type_cits == 3 and avg_point == 1:
            rho_pre1 = calc_rho(rho, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pre2 = calc_rho(rho, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pp1 = calc_rho(rho, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)

            ct1 = (pp_diffkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * rho_pre1 / (ssqtkk ** 2)))
            ct2 = 2 * pp_diffkk * bkk * rho_pre2 / ssqtkk
            ct3 = 2 * pp_diffkk * bkk * rho_pp1 / ssqtkk
            cits_theta = icc * (ct1 + ct2 - ct3)

            ce4 = (pp_diffkk ** 2) / ssqtkk
            cits_e = ((1 - icc) / nsamp) * ce4

            if did_cits == 2:
                cits_term1 += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
            its_term1 += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)

            if cross_sec == 2:
                phi_pre1 = calc_rho(phi, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pre2 = calc_rho(phi, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pp1 = calc_rho(phi, 8, 1, bkk, skk, P, tbarkk, 0, ar1, tpp)

                ce1 = (pp_diffkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * phi_pre1 / (ssqtkk ** 2)))
                ce2 = 2 * pp_diffkk * bkk * phi_pre2 / ssqtkk
                ce3 = 2 * pp_diffkk * bkk * phi_pp1 / ssqtkk
                cits_e_long = ((1 - icc) / nsamp) * (ce1 + ce2 - ce3)

                if did_cits == 2:
                    cits_term1_long += (akk ** 2) * ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                its_term1_long += (akk ** 2) * (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)

        # ---- CITS discrete, specific post ----
        if did_cits > 1 and type_cits == 3 and avg_point == 2 and iqkk == 1:
            rho_pre1 = calc_rho(rho, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pre2 = calc_rho(rho, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
            rho_pp1 = calc_rho(rho, 8, 1, bkk, qkk, qkk, tbarkk, 0, ar1, tpp)

            diff_qkk = tpp[qkk] - tbarkk

            ct1 = (diff_qkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * rho_pre1 / (ssqtkk ** 2)))
            ct2 = 2 * diff_qkk * bkk * rho_pre2 / ssqtkk
            ct3 = 2 * diff_qkk * bkk * rho_pp1 / ssqtkk
            cits_theta = icc * (ct1 + ct2 - ct3)

            ce4 = (diff_qkk ** 2) / ssqtkk
            cits_e = ((1 - icc) / nsamp) * ce4

            if did_cits == 2:
                cits_term1 += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e)
            its_term1 += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e)

            if cross_sec == 2:
                phi_pre1 = calc_rho(phi, 6, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pre2 = calc_rho(phi, 7, 1, bkk, 0, 0, tbarkk, 0, ar1, tpp)
                phi_pp1 = calc_rho(phi, 8, 1, bkk, qkk, qkk, tbarkk, 0, ar1, tpp)

                ce1 = (diff_qkk ** 2) * ((1 / ssqtkk) + ((bkk - 1) * bkk * phi_pre1 / (ssqtkk ** 2)))
                ce2 = 2 * diff_qkk * bkk * phi_pre2 / ssqtkk
                ce3 = 2 * diff_qkk * bkk * phi_pp1 / ssqtkk
                cits_e_long = ((1 - icc) / nsamp) * (ce1 + ce2 - ce3)

                if did_cits == 2:
                    cits_term1_long += ((1.0 / mtkk) + (1.0 / mckk)) * r2 * (cits_theta + cits_e_long)
                its_term1_long += (1.0 / mtkk_its) * r2 * (cits_theta + cits_e_long)

    # ---- AGGREGATE TOTAL VARIANCES ----
    did_tot = None
    cits_tot = None
    its_tot = None

    sumak2 = sumwk ** 2
    sumiq2 = sumiq ** 2 if sumiq > 0 else 1

    # DID total
    if cross_sec == 1:
        if avg_point == 1:
            did_tot = deff_wgt * (1.0 / sumak2) * term1
        else:
            did_tot = deff_wgt * (1.0 / sumiq2) * term1q
    else:
        if avg_point == 1:
            did_tot = deff_wgt * (1.0 / sumak2) * term1_long
        else:
            did_tot = deff_wgt * (1.0 / sumiq2) * term1q_long

    # CITS/ITS totals
    if did_cits > 1 and type_cits == 1:
        # Fully interacted: adds DID + CITS terms
        if cross_sec == 1:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * (term1 + cits_term1)
                its_tot = deff_wgt * (1.0 / sumak2) * (term1_noc + its_term1)
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * (term1 + cits_term1)
                its_tot = deff_wgt * (1.0 / sumiq2) * (term1_noc + its_term1)
        else:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * (term1_long + cits_term1_long)
                its_tot = deff_wgt * (1.0 / sumak2) * (term1_long_noc + its_term1_long)
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * (term1_long + cits_term1_long)
                its_tot = deff_wgt * (1.0 / sumiq2) * (term1_long_noc + its_term1_long)

    elif did_cits > 1 and type_cits == 2:
        # Common slopes: only CITS terms (no DID added)
        if cross_sec == 1:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * cits_term1
                its_tot = deff_wgt * (1.0 / sumak2) * its_term1
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * cits_term1
                its_tot = deff_wgt * (1.0 / sumiq2) * its_term1
        else:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * cits_term1_long
                its_tot = deff_wgt * (1.0 / sumak2) * its_term1_long
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * cits_term1_long
                its_tot = deff_wgt * (1.0 / sumiq2) * its_term1_long

    elif did_cits > 1 and type_cits == 3:
        # Discrete: DID + CITS terms
        if cross_sec == 1:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * (term1 + cits_term1)
                its_tot = deff_wgt * (1.0 / sumak2) * (term1_noc + its_term1)
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * (term1q + cits_term1)
                its_tot = deff_wgt * (1.0 / sumiq2) * (term1q_noc + its_term1)
        else:
            if avg_point == 1:
                cits_tot = deff_wgt * (1.0 / sumak2) * (term1_long + cits_term1_long)
                its_tot = deff_wgt * (1.0 / sumak2) * (term1_long_noc + its_term1_long)
            else:
                cits_tot = deff_wgt * (1.0 / sumiq2) * (term1q_long + cits_term1_long)
                its_tot = deff_wgt * (1.0 / sumiq2) * (term1q_long_noc + its_term1_long)

    return {
        'did_tot': did_tot,
        'cits_tot': cits_tot,
        'its_tot': its_tot,
        'sumak': sumak,
        'sumiq': sumiq,
    }


def get_design_label(inputs):
    """Return a human-readable label for the design configuration."""
    d = inputs.design_type
    c = inputs.cits_spec
    p = inputs.panel_type

    design_names = {
        (DesignType.DID, None, PanelType.CROSS_SECTIONAL): "DID Cross-Sectional Design",
        (DesignType.DID, None, PanelType.LONGITUDINAL): "DID Longitudinal Design",
    }

    if d == DesignType.DID:
        panel = "Cross-Sectional" if p == PanelType.CROSS_SECTIONAL else "Longitudinal"
        return f"DID {panel} Design"
    elif d == DesignType.CITS:
        panel = "Cross-Sectional" if p == PanelType.CROSS_SECTIONAL else "Longitudinal"
        spec = {1: "Fully-Interacted", 2: "Common-Slopes", 3: "Discrete"}[int(c)]
        return f"CITS {spec} {panel} Design"
    elif d == DesignType.ITS:
        panel = "Cross-Sectional" if p == PanelType.CROSS_SECTIONAL else "Longitudinal"
        spec = {1: "Fully-Interacted", 2: "Common-Slopes", 3: "Discrete"}[int(c)]
        return f"ITS {spec} {panel} Design"
    return "Unknown Design"


def run_power_analysis(inputs):
    """
    Main entry point. Returns results dict.
    """
    time_info = setup_time_periods(inputs)
    var_result = compute_variance(inputs, time_info)

    did_cits = int(inputs.design_type)
    type_cits = int(inputs.cits_spec)
    avg_point = int(inputs.post_period_type)

    # Select the appropriate variance
    if did_cits == 1:
        var_tot = var_result['did_tot']
    elif did_cits == 2:
        var_tot = var_result['cits_tot']
    elif did_cits == 3:
        var_tot = var_result['its_tot']

    if var_tot is None or var_tot <= 0:
        return {
            'results': [],
            'design_label': get_design_label(inputs),
            'variance': var_tot,
            'error': 'Computed variance is non-positive. Check inputs.'
        }

    # Determine total clusters for DF
    if inputs.analysis_mode == AnalysisMode.CALC_MDE:
        if did_cits < 3:
            mclus = inputs.mt + inputs.mc
        else:
            mclus = inputs.mt

    # Power range
    results = []
    power_vals = np.arange(inputs.power_min, inputs.power_max + inputs.power_step / 2, inputs.power_step)

    sumak = var_result['sumak']
    sumiq = var_result['sumiq']
    P = inputs.n_time_periods
    K = inputs.n_timing_groups

    if inputs.analysis_mode == AnalysisMode.CALC_MDE:
        for pwr in power_vals:
            if avg_point == 2:
                m_for_df = mclus  # approximate
            else:
                m_for_df = mclus

            df = compute_df(inputs.design_type, inputs.cits_spec,
                            m_for_df, P, K, sumak, sumiq, avg_point)
            mde_val = compute_mde(var_tot, df, inputs.alpha, pwr, inputs.two_tailed)
            results.append({'power': round(pwr, 2), 'mde_value': round(mde_val, 4)})

    else:  # CALC_CLUSTERS
        for pwr in power_vals:
            def df_func(m):
                return compute_df(inputs.design_type, inputs.cits_spec,
                                  m, P, K, sumak, sumiq, avg_point)

            opt = compute_required_clusters(
                var_tot, inputs.target_mde, inputs.alpha, pwr,
                inputs.two_tailed, df_func
            )
            if opt['converged']:
                results.append({'power': round(pwr, 2), 'required_clusters': opt['m_opt']})
            else:
                results.append({'power': round(pwr, 2), 'required_clusters': None})

    return {
        'results': results,
        'design_label': get_design_label(inputs),
        'variance': var_tot,
    }
