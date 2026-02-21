"""
Scenario comparison: generate MDE grids across ICC and R-squared values.
"""

import copy
import numpy as np
import pandas as pd
from config import PowerInputs, AnalysisMode
from power_engine import run_power_analysis


def generate_scenario_grid(base_inputs, icc_values, r2yx_values, fixed_power=0.80):
    """
    Compute MDE (or required clusters) for each ICC x R2yx combination.

    Returns a pandas DataFrame with ICC as rows, R2yx as columns.
    """
    n_icc = len(icc_values)
    n_r2 = len(r2yx_values)
    results = np.full((n_icc, n_r2), np.nan)

    is_mde_mode = base_inputs.analysis_mode == AnalysisMode.CALC_MDE

    for i, icc_val in enumerate(icc_values):
        for j, r2_val in enumerate(r2yx_values):
            modified = copy.deepcopy(base_inputs)
            modified.icc = icc_val
            modified.r2yx = r2_val
            modified.power_min = fixed_power
            modified.power_max = fixed_power
            modified.power_step = 0.05

            result = run_power_analysis(modified)

            if result['results']:
                row = result['results'][0]
                if is_mde_mode:
                    results[i, j] = row.get('mde_value', np.nan)
                else:
                    val = row.get('required_clusters')
                    results[i, j] = val if val is not None else np.nan

    row_labels = [f"{v:.2f}" for v in icc_values]
    col_labels = [f"{v:.2f}" for v in r2yx_values]

    df = pd.DataFrame(results, index=row_labels, columns=col_labels)
    df.index.name = "ICC"
    df.columns.name = "R2yx"

    return df
