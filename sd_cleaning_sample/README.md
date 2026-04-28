# sd_cleaning (sample: baseline)

**Minimal** code that isolates **one step** of the snow-depth pipeline: baseline correction in low-variance segments (`detect_flat_offsets` in `sd_ground_baseline.py`). The rest of the processing lives outside this folder.

## Contents

| File | Role |
|------|------|
| `sd_ground_baseline.py` | Algorithm + optional CLI for visualization |
| `baseline_viz.ipynb` | Same logic, step by step, in Jupyter |
| `SD_ValleOlivares_05706003.csv` | Reproducible example (`date` + `05706003` series) |
| `requirements.txt` | Dependencies |

## Quick start

```bash
python -m pip install -r requirements.txt
```

- **Notebook**: open `baseline_viz.ipynb` and run the cells in order (adjust the CSV path if needed).
- **Command line** with the default sample CSV:

```bash
python sd_ground_baseline.py
```

Other CSVs need a `date` column and the values column (the only numeric column besides date, or use `--col`). Without a graphical environment, you can save the figure with `-o output.png`.
