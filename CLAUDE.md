# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A Streamlit web app for statistical power analysis of panel data designs, translated from the R Shiny app in Schochet (2022, JEBS). Supports DID, CITS, and ITS designs with multiple specifications.

## Commands

```bash
# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
# Headless mode (e.g., remote server)
streamlit run app.py --server.headless true

# Run verification benchmarks (no test framework — standalone script)
python verify_power.py
```

There is no linter, formatter, or test framework configured.

## Architecture

The app follows a strict separation: **`power_engine.py` has zero Streamlit imports** and is pure Python/NumPy/SciPy. All UI lives in `app.py`.

### Module Responsibilities

- **`config.py`** — Enums (`DesignType`, `CITSSpec`, `AnalysisMode`, `PanelType`, `AutocorrStructure`, etc.) and the `PowerInputs` dataclass that carries all parameters to the engine.
- **`validation.py`** — `validate_inputs(PowerInputs) -> ValidationResult` with error messages.
- **`power_engine.py`** (~1000 lines) — Core computation: `run_power_analysis(inputs)` is the entry point. Internally calls `setup_time_periods` → `compute_variance` → `compute_mde` or `compute_required_clusters` per power level.
- **`scenarios.py`** — `generate_scenario_grid()` sweeps ICC × R² combinations calling `run_power_analysis` for each cell.
- **`app.py`** (~530 lines) — Streamlit UI with sidebar inputs (5 expander sections) and 3 main tabs (Results, Scenario Comparison, Help).
- **`verify_power.py`** — Benchmarks against Schochet (2022) Table 3 and independent formula implementations.

### Data Flow

```
Sidebar widgets → PowerInputs dataclass → validate_inputs() → run_power_analysis()
  → setup_time_periods() → compute_variance() → compute_mde/compute_required_clusters per power level
  → DataFrame → Plotly chart + table + CSV download
```

## Critical Implementation Details

- **1-based indexing convention**: The R original uses 1-based arrays. The Python translation uses **padded NumPy arrays where index 0 is unused**, faithfully mirroring R indexing. Do not "fix" this to 0-based — it would break alignment with the R source and verification benchmarks.
- **Variable naming from R**: The R variable `r2` (around line 1514 of the `.r` file) is actually `(1-r2yx)/(1-r2tx)` — a precision factor, not an R-squared. Variable names follow the R source for traceability.
- **Reference R file**: `JEBS_Power_Panel_Dashboard_Schochet.r` (~3150 lines) is the original R Shiny app kept for reference.
- **Verification**: All 6 benchmarks from Schochet (2022) Table 3 match exactly; independent formula checks match within < 0.0001.
