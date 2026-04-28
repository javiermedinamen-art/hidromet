<p align="center">
  <img src="snowflake-favicon.svg" width="96" height="96" alt="Snow Depth Observatory — snowflake mark">
</p>

<h1 align="center">Snow Depth Observatory</h1>

<p align="center">
  <strong>Southern Andes</strong> · curated snow-depth observations, stations, and watershed context
</p>

<p align="center">
  <em>A quiet, map-first workspace to explore quality-controlled daily records turned into monthly and seasonal metrics—without losing sight of where each measurement lives in the landscape.</em>
</p>

---

## What this is

This repository (**HidroMet**) is a **static web front end** for **Snow Depth Observatory**: a small scientific portal aimed at researchers, water managers, and anyone who needs a clear view of **snow depth** across the **Southern Andes**. It connects **station metadata and time series** to **basin and sub-basin boundaries** so you can move from “what does this site show?” to “how does this watershed behave?” in one place.

The live experience is intentionally **readable and visual**: an interactive map (Leaflet), **Plotly** charts for station series and **watershed comparisons**, and short methodological notes so the aggregates stay interpretable.

## What you can do here

- **Browse stations** on the map, with filters for basin hierarchy and intuitive search.
- **Open a station** to inspect its snow-depth time series and related diagnostics.
- **Compare watersheds** by selecting a basin or sub-basin polygon—monthly and **annual (May–October)** views summarize conditions under documented coverage rules.
- **Follow the methodology** embedded in the UI (e.g. minimum daily coverage for monthly inclusion, rules for annual medians).

## Repository layout (high level)

| Piece | Role |
|--------|------|
| `index.html` | Main app: map, charts, onboarding, and assets wired for deployment |
| `export_data.py` | Python pipeline: builds the static CSV/JSON/GeoJSON consumed by the app |
| `data_static/` | Generated artifacts used by the browser (after running the export script) |
| `snowflake-favicon.svg` | Icon / brand mark used across the site |

## Running locally

Serve the project root over HTTP (many features expect fetches rather than `file://`). For example:

```bash
python -m http.server 8080
```

Then open `http://localhost:8080/` and ensure `data_static` (or equivalent outputs) is present if the map requires those assets.

---

<p align="center">
  Built for transparent snow hydrology in the mountains between sky and catchment.
</p>
