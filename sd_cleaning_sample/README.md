# Southern Andes snow-depth QC — reproducibility sample

This folder is a **reduced reproducibility package** for the quality-control (QC) workflow described in the manuscript. It does **not** redistribute the full operational pipeline.

Analyses in the paper used **81 stations**. The open daily dataset redistributes **69 stations** because of institutional restrictions on some series (e.g. UdeChile). Parameter tables below still document the basin-level settings used when producing the full QC’d product.

## What is included vs not included

| Manuscript step | In this package? |
|-----------------|------------------|
| Step 1 — Compilation of multi-source SD | No (described in the manuscript) |
| **Step 2.1 — Zero-level / baseline correction** | **Yes** (`detect_flat_offsets`) |
| Step 2.2 — Spike removal (MAD) | No (parameters documented in the table) |
| Step 2.3 — Physical range check | No (parameters documented in the table) |
| Step 3 — Harmonization of overlapping series | No (described in the manuscript) |
| Step 4 — Cross-variable validation (PCI) | No (described in the manuscript) |

**Automated processing order** (Step 2):  
`detect_flat_offsets` → `detect_spikes_mad` → `apply_physical_range`, with optional station-level expert filters before/after. Only Step 2.1 is executable here.

## Contents

| File | Role |
|------|------|
| `sd_ground_baseline.py` | Step 2.1 implementation + CLI |
| `baseline_viz.ipynb` | Same logic, interactive |
| `SD_ValleOlivares_05706003.csv` | Example input (Valle Olivares, Maipo, 33° S) |
| `QC_parameters_by_river_basin.csv` | Basin-level parameters for Steps 2.1–2.3 |
| `requirements.txt` | Python dependencies |

## Expected input / output formats

**Input (sample):** CSV with

- `date` — daily timestamps (`YYYY-MM-DD`)
- one numeric snow-depth column in **cm** (here `05706003`)

**Output of `detect_flat_offsets`:**

- corrected DataFrame (same shape; SD with piecewise baseline removed)
- offset series (cm) applied to each column

The published multi-station product keeps the same tabular layout as the raw compilation (`date` + one column per station), with quality-controlled values.

## Default parameters in this sample (Maipo / Valle Olivares)

These match the Maipo (DGA) block used in the operational notebook:

| Parameter | Value | Role |
|-----------|-------|------|
| `window` | 5 days | Moving max−min window for stable segments |
| `var_thr` | 0.5 cm | Max range to call a segment stable |
| `tol` | 0.3 cm | Near-zero offsets set to 0 |
| `drift_thr` | 0.3 cm | Minimum change to update the baseline |
| `snow_thr` | 5.0 cm | Do not update baseline if median looks like snow |

Basin-specific values for **all** river basins (including spike and physical-range settings, with `min_cm = 0` everywhere) are in `QC_parameters_by_river_basin.csv`.

## Quick start

```bash
python -m pip install -r requirements.txt
```

**Notebook:** open `baseline_viz.ipynb` and run cells in order.

**CLI** (default sample CSV):

```bash
python sd_ground_baseline.py
```

Save a figure without a display:

```bash
python sd_ground_baseline.py -o baseline_valle_olivares.png
```

Other single-station CSVs need a `date` column and one numeric SD column (or pass `--col`).

## Suggested Code availability text (manuscript §5)

```text
The full quality-control workflow comprises four sequential steps described in
Section 2 (compilation; automated preprocessing with zero-level correction,
spike removal, and physical-range filtering; series harmonization; and
cross-variable validation).

A reduced, executable sample of Step 2.1 (zero-level / ground-reference
correction; detect_flat_offsets) is available at
https://github.com/javiermedinamen-art/hidromet/tree/main/sd_cleaning_sample,
applied to the Valle Olivares station (05706003; 33° S). The repository
includes the Python module, a Jupyter notebook, an example input CSV,
dependency specifications, and a basin-level parameter table for Steps 2.1–2.3
(QC_parameters_by_river_basin). Steps 2.2, 2.3, 3, and 4 are fully specified in
Section 2 of the manuscript (algorithms, decision rules, and thresholds) but
are not redistributed as source code. The open dataset provides quality-
controlled daily snow depth for 69 stations; institutional restrictions prevent
redistribution of the remaining series used in the 81-station analysis.
```
