"""
Configuration module: Enums, dataclass, and defaults for the Power Panel app.
Based on Schochet (JEBS, 2022) power analysis for panel data designs.
"""

from enum import IntEnum
from dataclasses import dataclass, field
from typing import List, Optional


class DesignType(IntEnum):
    DID = 1
    CITS = 2
    ITS = 3


class CITSSpec(IntEnum):
    FULLY_INTERACTED = 1   # Different pre- and post-period slopes
    COMMON_SLOPES = 2      # Common pre- and post-period slopes
    DISCRETE = 3           # Discrete post-period indicators (like DID)


class AnalysisMode(IntEnum):
    CALC_MDE = 1           # Calculate MDE given sample size
    CALC_CLUSTERS = 2      # Calculate required clusters given target MDE


class PanelType(IntEnum):
    CROSS_SECTIONAL = 1
    LONGITUDINAL = 2


class AutocorrStructure(IntEnum):
    AR1 = 1
    CONSTANT = 2
    NONE = 3


class PostPeriodType(IntEnum):
    AVERAGE = 1
    SPECIFIC = 2


class SpecificPostUnit(IntEnum):
    TIME_PERIOD = 1
    EXPOSURE_POINT = 2


@dataclass
class PowerInputs:
    design_type: DesignType = DesignType.DID
    cits_spec: CITSSpec = CITSSpec.FULLY_INTERACTED
    analysis_mode: AnalysisMode = AnalysisMode.CALC_MDE
    panel_type: PanelType = PanelType.CROSS_SECTIONAL

    n_time_periods: int = 10
    evenly_spaced: bool = True
    time_intervals: Optional[List[int]] = None

    n_timing_groups: int = 2
    start_times: List[int] = field(default_factory=lambda: [4, 6])

    post_period_type: PostPeriodType = PostPeriodType.AVERAGE
    specific_post_unit: SpecificPostUnit = SpecificPostUnit.TIME_PERIOD
    q_time: int = 6
    l_time: int = 6

    # Sample size inputs (mode 1: given sample)
    mt: int = 20
    mtk: List[int] = field(default_factory=lambda: [10, 10])
    mc: int = 20
    mck: List[int] = field(default_factory=lambda: [10, 10])

    # Sample size inputs (mode 2: find clusters)
    target_mde: float = 0.20
    rt: float = 0.50
    rtk: List[float] = field(default_factory=lambda: [0.50, 0.50])
    rck: List[float] = field(default_factory=lambda: [0.50, 0.50])

    n_per_cluster: int = 100

    # Error structure
    autocorr_structure: AutocorrStructure = AutocorrStructure.AR1
    rho: float = 0.40
    phi: float = 0.40
    icc: float = 0.05

    # Precision gains
    r2yx: float = 0.0
    r2tx: float = 0.0
    deff_wgt: float = 1.0

    # Hypothesis testing
    alpha: float = 0.05
    two_tailed: bool = True
    power_min: float = 0.60
    power_max: float = 0.90
    power_step: float = 0.05
