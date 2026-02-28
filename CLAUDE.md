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
- **`validation.py`** — `validate_inputs(PowerInputs) -> ValidationResult`. `ValidationResult` has `is_valid: bool`, `errors: List[str]`, and `add_error(msg)`.
- **`power_engine.py`** (~1000 lines) — Core computation: `run_power_analysis(inputs)` is the entry point. Internally calls `setup_time_periods` → `compute_variance` → `compute_mde` or `compute_required_clusters` per power level. Key internal functions:
  - `calc_rho()` — autocorrelation calculator with 12 cases + 3 special cases (90, 100, 110); AR1 uses `rho^time_distance`, constant uses flat `rho`
  - `compute_df()` — degrees of freedom vary by design/spec/mode (e.g. DID avg: `M*P - M - K*P - sum_ak`; CITS fully-interacted avg: `M*P - 8*K`)
  - `compute_required_clusters()` — secant method, max 25 iterations, tol 1e-6; returns `{'m_opt': int, 'converged': bool}` or `None`
- **`scenarios.py`** — `generate_scenario_grid(base_inputs, icc_values, r2yx_values, fixed_power=0.80) -> pd.DataFrame`. Grid rows = ICC, columns = R2yx, cells = MDE or clusters at fixed power.
- **`app.py`** (~640 lines) — Streamlit UI with sidebar inputs (5 expander sections) and 3 main tabs (Results, Scenario Comparison, Help). The Help tab has 8 expanders, including a "Statistical Formulas" section with rendered LaTeX for all variance, MDE, and df formulas.
- **`verify_power.py`** — Benchmarks against Schochet (2022) Table 3 and independent formula implementations.

### Data Flow

```
Sidebar widgets → PowerInputs dataclass → validate_inputs() → run_power_analysis()
  → setup_time_periods() → compute_variance() → compute_mde/compute_required_clusters per power level
  → DataFrame → Plotly chart + table + CSV download
```

### Engine Return Shapes (all functions return dicts)

- `setup_time_periods()` → `{tpp, tbar, ssqt, tbara, ssqta, pp_diff, bk, ak, tbar_full, ssqt_full}`
- `compute_variance()` → `{did_tot, cits_tot, its_tot, sumak, sumiq}`
- `run_power_analysis()` → `{results, design_label, variance, error?}`

## Critical Implementation Details

- **1-based indexing convention**: The R original uses 1-based arrays. The Python translation uses **padded NumPy arrays where index 0 is unused**, faithfully mirroring R indexing. Do not "fix" this to 0-based — it would break alignment with the R source and verification benchmarks.
- **Variable naming from R**: The R variable `r2` (around line 1514 of the `.r` file) is actually `(1-r2yx)/(1-r2tx)` — a precision factor, not an R-squared. Variable names follow the R source for traceability.
- **Reference R file**: `JEBS_Power_Panel_Dashboard_Schochet.r` (~3150 lines) is the original R Shiny app kept for reference.
- **Verification**: All 6 benchmarks from Schochet (2022) Table 3 match exactly; independent formula checks match within < 0.0001.
- **Streamlit LaTeX rendering**: `st.markdown()` supports `$...$` (inline) and `$$...$$` (display) math. Do not put a newline inside a `$$` block — Streamlit closes the block early and the remainder renders as raw LaTeX text. Each display formula must be a single unbroken line.
