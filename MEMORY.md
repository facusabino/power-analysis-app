# Power Panel Project

## What it is
A Streamlit Python app for power analysis of panel data designs (DID, CITS, ITS), based on Schochet (JEBS, 2022). Replaces an existing R Shiny app that was clunky.

## Key files
- `config.py` — Enums + `PowerInputs` dataclass
- `power_engine.py` — Core computation (calc_rho, compute_variance, compute_mde, etc.)
- `validation.py` — Input validation
- `scenarios.py` — ICC x R2 scenario grid generation
- `app.py` — Streamlit UI (sidebar + Results/Scenario tabs)
- `verify_power.py` — Independent verification script (benchmarks against paper)
- `requirements.txt` — streamlit, numpy, scipy, pandas, plotly

## Reference files
- `JEBS_Power_Panel_Dashboard_Schochet.r` — Original R Shiny app (~3150 lines)
- `Schochet paper.pdf` — Academic paper with formulas
- `power panel app.pdf` — App documentation

## Architecture
- `power_engine.py` has zero Streamlit dependency (pure computation, testable)
- R code uses 1-based indexing; Python uses padded arrays (index 0 = unused) for faithful translation
- R's `r2` variable (line 1514) is actually `(1-r2yx)/(1-r2tx)` — a precision factor, not R-squared

## Verification
- All 6 benchmarks from Schochet (2022) Table 3 match exactly
- Independent formula implementation matches engine within < 0.0001 across 8 test cases
- Run `python verify_power.py` to re-run the audit

## Run command
```
streamlit run app.py --server.headless true
```
Opens at http://localhost:8501
