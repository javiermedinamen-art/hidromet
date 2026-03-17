import json
import os
from datetime import datetime
from difflib import SequenceMatcher
import unicodedata

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

OUTPUT_DIR = "data_static"
PATH_ESTACIONES = "data/estaciones.csv"
SD_FILE_PATH = "data/SD_clean_v3.csv"
SD_UNIT = "cm"
SD_DISPLAY_NAME = "Snow_Depth"
FILL_RULE_PERCENT = 0.7
MONTHLY_AGG_METHOD = "median"
MONTHS_REQUIRED = [5, 6, 7, 8, 9, 10]

PATH_CUENCAS_SHP = "data/cuencas_bna/cuencas_bna.shp"
ID_COLUMNA_CUENCA = "COD_CUEN"
PATH_CUENCAS_ARG_SHP = "data/cuencas_arg/cuencas_arg.shp"
PATH_SUBCUENCAS_SHP = "data/subcuencas_bna/subcuencas_bna.shp"
ID_COLUMNA_SUBCUENCA = "COD_SUBC"


def cleanup_output_dir():
    """
    Remove all previously generated data files so each run produces a full refresh.
    Data in 'data/' (estaciones, SD series, shapefiles) changes over time;
    this ensures data_static always reflects the current sources.
    """
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        return

    removed = 0
    for filename in os.listdir(OUTPUT_DIR):
        path = os.path.join(OUTPUT_DIR, filename)
        if os.path.isfile(path) and filename.endswith((".csv", ".json", ".geojson")):
            os.remove(path)
            removed += 1
    if removed:
        print(f"Cleaned {removed} previous output file(s) from {OUTPUT_DIR}.")


def cargar_hierarquia(path):
    """Read a hierarchy shapefile and return it in EPSG:4326."""
    if not os.path.exists(path):
        print(f"Warning: {path} not found. Skipping.")
        return None
    try:
        return gpd.read_file(path, encoding="utf-8").to_crs("EPSG:4326")
    except Exception:
        try:
            return gpd.read_file(path, encoding="latin-1").to_crs("EPSG:4326")
        except Exception as exc:
            print(f"Critical error reading {path}: {exc}")
            return None
    except Exception as exc:
        print(f"Critical error reading {path}: {exc}")
        return None


def cargar_cuencas_chile_estandar():
    """Load Chilean basins with standard basin ID/name columns."""
    gdf = cargar_hierarquia(PATH_CUENCAS_SHP)
    if gdf is None or gdf.empty:
        return gpd.GeoDataFrame(columns=["COD_CUEN", "NOM_CUEN", "geometry"], crs="EPSG:4326")
    return gdf[["COD_CUEN", "NOM_CUEN", "geometry"]].copy()


def cargar_cuencas_argentina_estandar():
    """Load Argentine basins and map them to the same basin schema used by the app."""
    gdf = cargar_hierarquia(PATH_CUENCAS_ARG_SHP)
    if gdf is None or gdf.empty:
        return gpd.GeoDataFrame(columns=["COD_CUEN", "NOM_CUEN", "geometry"], crs="EPSG:4326")

    gdf = gdf[["CUENCA", "geometry"]].copy()
    gdf = gdf[gdf["CUENCA"].notna()].copy()
    gdf["CUENCA"] = gdf["CUENCA"].astype(str).str.strip()
    gdf = gdf[gdf["CUENCA"] != ""].copy()
    gdf = gdf.drop_duplicates(subset=["CUENCA"]).reset_index(drop=True)
    gdf["COD_CUEN"] = [f"ARG_{idx:03d}" for idx in range(1, len(gdf) + 1)]
    gdf["NOM_CUEN"] = gdf["CUENCA"]
    return gdf[["COD_CUEN", "NOM_CUEN", "geometry"]].copy()


def obtener_cuencas_combinadas():
    """Combine Chilean and Argentine basins into one basin layer."""
    chile = cargar_cuencas_chile_estandar()
    argentina = cargar_cuencas_argentina_estandar()
    frames = [gdf for gdf in (chile, argentina) if gdf is not None and not gdf.empty]
    if not frames:
        return gpd.GeoDataFrame(columns=["COD_CUEN", "NOM_CUEN", "geometry"], crs="EPSG:4326")
    return gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs="EPSG:4326")


def cargar_estaciones():
    """Read station metadata and return a GeoDataFrame."""
    if not os.path.exists(PATH_ESTACIONES):
        print(f"Critical error: stations file not found at {PATH_ESTACIONES}")
        return None

    df = None
    for encoding in ("utf-8", "cp1252", "latin-1"):
        try:
            df = pd.read_csv(
                PATH_ESTACIONES,
                encoding=encoding,
                sep=",",
                decimal=".",
                dtype={"code_internal": str},
            )
            print(f"Reading stations from {PATH_ESTACIONES} with encoding {encoding}")
            break
        except Exception:
            continue

    if df is None:
        print(f"Critical error: could not read {PATH_ESTACIONES}")
        return None

    df.columns = [c.strip().lower() for c in df.columns]
    lon_col = "lon" if "lon" in df.columns else "long" if "long" in df.columns else None
    required_cols = ["lat", "code_internal", "name", "fuente", "basin"]
    if lon_col is None or not all(col in df.columns for col in required_cols):
        print("Error: stations.csv is missing required columns.")
        print(f"Columns found: {list(df.columns)}")
        return None

    if "elevation" not in df.columns:
        df["elevation"] = np.nan

    df["longitude"] = pd.to_numeric(df[lon_col].astype(str).str.strip(), errors="coerce")
    df["latitude"] = pd.to_numeric(df["lat"].astype(str).str.strip(), errors="coerce")
    df["elevation"] = pd.to_numeric(df["elevation"], errors="coerce")
    df["code_internal"] = df["code_internal"].astype(str).str.strip().str.lower()
    df["name"] = df["name"].astype(str).str.strip()
    df["fuente"] = df["fuente"].where(df["fuente"].notna(), "").astype(str).str.strip()
    df["basin"] = df["basin"].where(df["basin"].notna(), "").astype(str).str.strip()

    df_clean = df.dropna(subset=["longitude", "latitude"]).copy()
    if df_clean.empty:
        print("Warning: no stations with valid coordinates remain.")
        return gpd.GeoDataFrame()

    geometry = [Point(xy) for xy in zip(df_clean["longitude"], df_clean["latitude"])]
    return gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")


def completar_cuencas_en_memoria(gdf_estaciones):
    """
    Fill missing basin names in memory using station coordinates and the
    combined Chile/Argentina basin layer, without editing estaciones.csv.
    """
    if gdf_estaciones is None or gdf_estaciones.empty:
        return gdf_estaciones

    gdf = gdf_estaciones.copy()
    basin_series = gdf["basin"].fillna("").astype(str).str.strip()
    missing_mask = basin_series.isin({"", "-", "nan", "NaN"})
    if not missing_mask.any():
        return gdf

    cuencas = obtener_cuencas_combinadas()
    if cuencas.empty:
        return gdf

    gdf_missing = gdf[missing_mask].copy()
    if gdf_missing.empty:
        return gdf

    joined = gpd.sjoin(gdf_missing, cuencas[["NOM_CUEN", "geometry"]], how="left", predicate="within")
    filled = 0
    for idx in gdf_missing.index:
        matches = joined.loc[joined.index == idx, "NOM_CUEN"].dropna().unique().tolist()
        if matches:
            gdf.at[idx, "basin"] = matches[0]
            filled += 1

    if filled:
        print(f"Filled {filled} missing basin name(s) in memory from basin shapefiles.")
    return gdf


def _normalize_code(code):
    """Normalize station code for matching: lowercase and remove spaces."""
    if pd.isna(code):
        return ""
    return str(code).strip().lower().replace(" ", "")


def _normalize_text_key(text):
    """Normalize free text for resilient station matching."""
    if pd.isna(text):
        return ""
    text = str(text).strip().lower()
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    return "".join(ch for ch in text if ch.isalnum())


def _codes_match(a, b):
    """Check if two codes match (handles leading zeros, spaces, case)."""
    na, nb = _normalize_code(a), _normalize_code(b)
    if na == nb:
        return True
    try:
        return int(na) == int(nb) if na.isdigit() and nb.isdigit() else False
    except (ValueError, TypeError):
        return False


def _text_match(a, b, threshold=0.86):
    """Fallback fuzzy match for text-like station identifiers."""
    na, nb = _normalize_text_key(a), _normalize_text_key(b)
    if not na or not nb:
        return False
    if na == nb:
        return True
    return SequenceMatcher(None, na, nb).ratio() >= threshold


def _name_first_word_match(est_name, sd_col):
    """Match when the first word of station name equals the SD column (case-insensitive)."""
    if pd.isna(est_name) or pd.isna(sd_col):
        return False
    first = str(est_name).strip().split()[0] if str(est_name).strip() else ""
    return _normalize_code(first) == _normalize_code(sd_col)


def build_station_mapping(gdf_estaciones, sd_columns):
    """
    Build mapping from SD column names (data IDs) to station metadata.
    Uses normalized matching (lowercase, no spaces, leading zeros) so
    'Piuquenes 10' matches 'Piuquenes10' and '430004' matches '0430004'.
    Also matches by first word of name (e.g. 'Portillo' matches 'Portillo Argentino').
    Returns a GeoDataFrame with SD column as code_internal.
    """
    estaciones = gdf_estaciones.copy()

    rows = []
    for sd_col in sd_columns:
        match = None
        for _, est_row in estaciones.iterrows():
            if _codes_match(est_row["code_internal"], sd_col):
                match = est_row
                break
        if match is None:
            for _, est_row in estaciones.iterrows():
                if _text_match(est_row["code_internal"], sd_col):
                    match = est_row
                    break
        if match is None:
            for _, est_row in estaciones.iterrows():
                if _name_first_word_match(est_row["name"], sd_col):
                    match = est_row
                    break
        if match is not None:
            row = match.copy()
            row["code_internal"] = sd_col
            rows.append(row)

    if not rows:
        return gpd.GeoDataFrame()

    return gpd.GeoDataFrame(rows, crs=estaciones.crs)


def cargar_series_temporales(file_path):
    """Load daily SD series and coerce all station columns to numeric."""
    print(f"\n--- Loading SD series from {file_path} ---")
    if not os.path.exists(file_path):
        print(f"Error: series file not found: {file_path}")
        return None

    df = None
    for encoding in ("utf-8", "cp1252", "latin-1"):
        try:
            df = pd.read_csv(file_path, encoding=encoding, sep=",", decimal=".")
            print(f"Reading series with encoding {encoding}")
            break
        except Exception:
            continue

    if df is None:
        print(f"Error: could not read {file_path}")
        return None

    df.columns = df.columns.astype(str).str.strip().str.lower()
    if "date" not in df.columns:
        print(f"Critical error: 'date' column not found in {file_path}")
        return None

    df["date"] = pd.to_datetime(df["date"].astype(str).str.strip(), errors="coerce")
    df.dropna(subset=["date"], inplace=True)
    df.set_index("date", inplace=True)

    for col in df.columns:
        df[col] = pd.to_numeric(df[col].astype(str).str.replace(",", ".", regex=False), errors="coerce")

    df.dropna(axis=1, how="all", inplace=True)
    df.dropna(axis=0, how="all", inplace=True)
    df.sort_index(inplace=True)
    print(f"Loaded SD series. Stations: {len(df.columns)}, rows: {len(df)}")
    return df


def calcular_agregacion_mensual(df_daily):
    """Compute valid monthly medians after applying the 70% daily coverage rule."""
    monthly_values = df_daily.resample("MS").agg(MONTHLY_AGG_METHOD)
    monthly_counts = df_daily.resample("MS").count()
    required_days = pd.Series(monthly_values.index.days_in_month * FILL_RULE_PERCENT, index=monthly_values.index)
    validity_mask = monthly_counts.ge(required_days, axis=0)
    monthly_values = monthly_values.where(validity_mask)
    monthly_values.dropna(axis=0, how="all", inplace=True)
    return monthly_values


def calcular_mediana_anual_mayo_octubre(df_daily, df_monthly_valid):
    """Compute annual May-Oct medians only for years with all required valid months."""
    monthly_required = df_monthly_valid[df_monthly_valid.index.month.isin(MONTHS_REQUIRED)]
    if monthly_required.empty:
        return pd.DataFrame(columns=df_daily.columns)

    valid_months_per_year = monthly_required.notna().groupby(monthly_required.index.year).sum()
    complete_years_mask = valid_months_per_year.eq(len(MONTHS_REQUIRED))

    daily_required = df_daily[df_daily.index.month.isin(MONTHS_REQUIRED)]
    annual_median = daily_required.groupby(daily_required.index.year).median()
    annual_median = annual_median.where(complete_years_mask.reindex(annual_median.index).fillna(False))
    annual_median.dropna(axis=0, how="all", inplace=True)
    return annual_median


def obtener_metadata_estaciones(gdf_estaciones):
    metadata = gdf_estaciones.drop_duplicates(subset=["code_internal"]).set_index("code_internal")
    return metadata


def normalizar_elevacion(valor):
    return None if pd.isna(valor) else float(valor)


def construir_payload_mensual(df_monthly, station_metadata):
    if df_monthly.empty:
        return {}

    complete_index = pd.date_range(df_monthly.index.min(), df_monthly.index.max(), freq="MS")
    df_complete = df_monthly.reindex(complete_index)
    payload = {}

    for station_code in df_complete.columns:
        if station_code not in station_metadata.index:
            continue
        station_info = station_metadata.loc[station_code]
        values = [None if pd.isna(val) else round(float(val), 1) for val in df_complete[station_code]]
        if any(value is not None for value in values):
            payload[station_code] = {
                "name": station_info["name"],
                "elevation": normalizar_elevacion(station_info.get("elevation")),
                "dates": df_complete.index.strftime("%Y-%m-%d").tolist(),
                "values": values,
            }

    return payload


def construir_payload_anual(df_annual, station_metadata):
    if df_annual.empty:
        return {}

    payload = {}
    years = [int(year) for year in df_annual.index.tolist()]
    for station_code in df_annual.columns:
        if station_code not in station_metadata.index:
            continue
        station_info = station_metadata.loc[station_code]
        values = [None if pd.isna(val) else round(float(val), 1) for val in df_annual[station_code]]
        if any(value is not None for value in values):
            payload[station_code] = {
                "name": station_info["name"],
                "elevation": normalizar_elevacion(station_info.get("elevation")),
                "years": years,
                "values": values,
            }

    return payload


def exportar_series_estacion(df_daily, df_annual, gdf_estaciones_sd):
    """Export daily and annual SD series for each station shown in the app."""
    print("\n--- Exporting station time series ---")
    valid_codes = set(gdf_estaciones_sd["code_internal"])
    exported_daily = 0
    exported_annual = 0

    for code in df_daily.columns:
        if code not in valid_codes:
            continue

        daily_series = pd.to_numeric(df_daily[code], errors="coerce").dropna()
        if not daily_series.empty:
            daily_output = pd.DataFrame({"Snow_Depth_cm": daily_series})
            daily_output.to_csv(
                os.path.join(OUTPUT_DIR, f"{code}_sd.csv"),
                index=True,
                index_label="Date",
                date_format="%Y-%m-%d",
                encoding="utf-8",
            )
            exported_daily += 1

        if code in df_annual.columns:
            annual_series = pd.to_numeric(df_annual[code], errors="coerce").dropna()
            if not annual_series.empty:
                annual_output = pd.DataFrame(
                    {
                        "Year": annual_series.index.astype(int),
                        "Snow_Depth_cm": annual_series.round(1).values,
                    }
                )
                annual_output.to_csv(
                    os.path.join(OUTPUT_DIR, f"{code}_sd_annual.csv"),
                    index=False,
                    encoding="utf-8",
                )
                exported_annual += 1

    print(f"Exported {exported_daily} daily station CSV files.")
    print(f"Exported {exported_annual} annual station CSV files.")


def procesar_jerarquia_geoespacial(gdf_hierarchy, id_columna, output_prefix, gdf_estaciones_sd, df_monthly, df_annual):
    """Export hierarchy GeoJSON plus monthly and annual SD JSON payloads."""
    print(f"\n--- Processing hierarchy: {output_prefix} ---")
    if gdf_hierarchy is None or gdf_hierarchy.empty:
        print(f"Warning: {output_prefix} hierarchy is empty. Skipping.")
        return
    gdf_hierarchy = gdf_hierarchy.to_crs(gdf_estaciones_sd.crs)

    gdf_join = gpd.sjoin(gdf_estaciones_sd, gdf_hierarchy, how="inner", predicate="within")
    if gdf_join.empty:
        print(f"Warning: no SD stations intersect {output_prefix}.")
        return

    ids_con_estaciones = gdf_join[id_columna].unique()
    gdf_filtrado = gdf_hierarchy[gdf_hierarchy[id_columna].isin(ids_con_estaciones)].copy()
    station_counts = gdf_join.groupby(id_columna).size().rename("n_estaciones")
    gdf_filtrado = gdf_filtrado.merge(station_counts, left_on=id_columna, right_index=True, how="left")
    gdf_filtrado["n_estaciones"] = gdf_filtrado["n_estaciones"].fillna(0).astype(int)

    geojson_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_sd.geojson")
    gdf_filtrado.to_file(geojson_path, driver="GeoJSON", encoding="utf-8")
    print(f"Exported hierarchy GeoJSON: {geojson_path}")

    station_metadata = obtener_metadata_estaciones(gdf_join)
    count_polygons = 0

    for poly_id in ids_con_estaciones:
        polygon_stations = gdf_join[gdf_join[id_columna] == poly_id].copy()
        polygon_stations.sort_values(by="elevation", ascending=False, na_position="last", inplace=True)
        station_codes = [code for code in polygon_stations["code_internal"].tolist() if code in df_monthly.columns]
        if not station_codes:
            continue

        monthly_polygon = df_monthly[station_codes].dropna(axis=0, how="all")
        monthly_json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_sd.json")
        annual_json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_sd_annual_median.json")

        monthly_payload = construir_payload_mensual(monthly_polygon, station_metadata)
        if monthly_payload:
            with open(monthly_json_path, "w", encoding="utf-8") as file_obj:
                json.dump(monthly_payload, file_obj, ensure_ascii=False, indent=2)

        annual_polygon = df_annual[station_codes].dropna(axis=0, how="all")
        annual_payload = construir_payload_anual(annual_polygon, station_metadata)
        if annual_payload:
            with open(annual_json_path, "w", encoding="utf-8") as file_obj:
                json.dump(annual_payload, file_obj, ensure_ascii=False, indent=2)

        if monthly_payload or annual_payload:
            count_polygons += 1

    print(f"Exported SD JSON files for {count_polygons} polygons in {output_prefix}.")


def exportar_datos_estaticos():
    print("Starting SnowData static export...")
    cleanup_output_dir()

    gdf_estaciones = cargar_estaciones()
    if gdf_estaciones is None or gdf_estaciones.empty:
        print("Station loading failed. Aborting.")
        return
    gdf_estaciones = completar_cuencas_en_memoria(gdf_estaciones)

    df_daily = cargar_series_temporales(SD_FILE_PATH)
    if df_daily is None or df_daily.empty:
        print("No valid SD series found. Aborting.")
        return

    station_codes_with_data = list(df_daily.columns)
    gdf_estaciones_sd = build_station_mapping(gdf_estaciones, station_codes_with_data)
    gdf_estaciones_sd["variables_disponibles"] = "SD"

    if gdf_estaciones_sd.empty:
        print("No station metadata matched the SD dataset. Aborting.")
        return

    mapped = set(gdf_estaciones_sd["code_internal"])
    codes_only_in_data = [c for c in station_codes_with_data if c not in mapped]
    if codes_only_in_data:
        print(f"Warning: {len(codes_only_in_data)} SD stations have no match in estaciones.csv")
        print(sorted(codes_only_in_data)[:10], "..." if len(codes_only_in_data) > 10 else "")

    df_monthly = calcular_agregacion_mensual(df_daily)
    df_annual = calcular_mediana_anual_mayo_octubre(df_daily, df_monthly)

    estaciones_geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    gdf_estaciones_sd[gdf_estaciones_sd.is_valid].to_file(
        estaciones_geojson_path, driver="GeoJSON", encoding="utf-8"
    )
    print(f"Exported stations GeoJSON: {estaciones_geojson_path}")

    station_names = (
        gdf_estaciones_sd.drop_duplicates(subset=["code_internal"])
        .set_index("code_internal")["name"]
        .to_dict()
    )
    with open(os.path.join(OUTPUT_DIR, "station_names.json"), "w", encoding="utf-8") as file_obj:
        json.dump(station_names, file_obj, ensure_ascii=False, indent=2)
    print("Exported station_names.json")

    exportar_series_estacion(df_daily, df_annual, gdf_estaciones_sd)

    cuencas_combinadas = obtener_cuencas_combinadas()
    subcuencas = cargar_hierarquia(PATH_SUBCUENCAS_SHP)

    procesar_jerarquia_geoespacial(cuencas_combinadas, ID_COLUMNA_CUENCA, "cuencas", gdf_estaciones_sd, df_monthly, df_annual)
    procesar_jerarquia_geoespacial(subcuencas, ID_COLUMNA_SUBCUENCA, "subcuencas", gdf_estaciones_sd, df_monthly, df_annual)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    version_for_cache = datetime.now().strftime("%Y-%m-%d-sd")
    source_order = ["DGA", "CEAZAMET", "IANIGLA", "UdeChile", "CIEP"]
    data_sources = sorted(
        [src for src in gdf_estaciones_sd["fuente"].dropna().unique().tolist() if src and src.lower() != "nan"],
        key=lambda source: source_order.index(source) if source in source_order else len(source_order),
    )

    metadata_path = os.path.join(OUTPUT_DIR, "export_metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "export_timestamp": timestamp,
                "data_version": version_for_cache,
                "data_sources": data_sources,
                "sources": {
                    "estaciones": PATH_ESTACIONES,
                    "sd_series": SD_FILE_PATH,
                },
            },
            f,
            indent=2,
        )
    print(f"Wrote {metadata_path} (data_version: {version_for_cache})")
    print("SnowData static export finished.")


if __name__ == "__main__":
    exportar_datos_estaticos()