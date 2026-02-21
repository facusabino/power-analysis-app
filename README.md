# Power Analysis App for Panel Data Designs

A Streamlit web application for conducting power analyses for commonly used panel data designs, based on the methods developed in [Schochet (2022)](https://doi.org/10.3102/10769986221112426).

## Supported Designs

- **Difference-in-Differences (DID)**
- **Comparative Interrupted Time Series (CITS)** — fully interacted, common slopes, and discrete specifications
- **Interrupted Time Series (ITS)** — fully interacted, common slopes, and discrete specifications

## Features

- **Two analysis modes**: calculate minimum detectable effect sizes (MDE) for a given sample, or calculate required cluster sample sizes for a target MDE
- **Scenario comparison**: generate tables of MDE values across grids of ICC and R-squared values — designed for easy copy-paste into reports
- **Flexible designs**: supports staggered treatment timing, AR(1) and constant autocorrelation structures, cross-sectional and longitudinal panels, uneven time spacing, and model covariates
- **Interactive UI**: collapsible parameter sections, dynamic timing group configuration, interactive Plotly charts, and CSV export

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

## Project Structure

| File | Description |
|------|-------------|
| `app.py` | Streamlit UI — sidebar inputs, results table/chart, scenario comparison tab |
| `power_engine.py` | Core computation module (no Streamlit dependency) |
| `config.py` | Enums and `PowerInputs` dataclass |
| `validation.py` | Input validation with structured error messages |
| `scenarios.py` | ICC x R-squared scenario grid generation |
| `verify_power.py` | Independent verification script — benchmarks against the paper |

## Verification

The computations have been independently verified against:

- **Table 3 of Schochet (2022)**: all 6 benchmark values for required cluster counts match exactly (DID, CITS fully-interacted, and CITS common-slopes across two parameter configurations)
- **Independent formula implementation**: MDE values from a second, standalone implementation of the paper's variance formulas match the engine within < 0.0001 across 8 test cases

Run the verification yourself:

```bash
python verify_power.py
```

## References

Schochet, P.Z. (2022). Statistical Power for Estimating Treatment Effects Using Difference-in-Differences and Comparative Interrupted Time Series Estimators with Variation in Treatment Timing. *Journal of Educational and Behavioral Statistics*, 48(6), 713-751.

## License

MIT
