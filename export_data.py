import json
import os
from datetime import datetime

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
PATH_SUBCUENCAS_SHP = "data/subcuencas_bna/subcuencas_bna.shp"
ID_COLUMNA_SUBCUENCA = "COD_SUBC"
PATH_SUBSUBCUENCAS_SHP = "data/subsubcuencas_bna/subsubcuencas_bna.shp"
ID_COLUMNA_SUBSUBCUENCA = "COD_SSUBC"


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


def cargar_estaciones():
    """Read station metadata and return a GeoDataFrame."""
    if not os.path.exists(PATH_ESTACIONES):
        print(f"Critical error: stations file not found at {PATH_ESTACIONES}")
        return None

    df = None
    for encoding in ("utf-8", "cp1252"):
        try:
            df = pd.read_csv(PATH_ESTACIONES, encoding=encoding, sep=",", decimal=".")
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
    df["fuente"] = df["fuente"].astype(str).str.strip()
    df["basin"] = df["basin"].astype(str).str.strip()

    df_clean = df.dropna(subset=["longitude", "latitude"]).copy()
    if df_clean.empty:
        print("Warning: no stations with valid coordinates remain.")
        return gpd.GeoDataFrame()

    geometry = [Point(xy) for xy in zip(df_clean["longitude"], df_clean["latitude"])]
    return gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")


def cargar_series_temporales(file_path):
    """Load daily SD series and coerce all station columns to numeric."""
    print(f"\n--- Loading SD series from {file_path} ---")
    if not os.path.exists(file_path):
        print(f"Error: series file not found: {file_path}")
        return None

    df = None
    for encoding in ("utf-8", "cp1252"):
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
                annual_output.to_csv(os.path.join(OUTPUT_DIR, f"{code}_sd_annual.csv"), index=False)
                exported_annual += 1

    print(f"Exported {exported_daily} daily station CSV files.")
    print(f"Exported {exported_annual} annual station CSV files.")


def procesar_jerarquia_geoespacial(shapefile_path, id_columna, output_prefix, gdf_estaciones_sd, df_monthly, df_annual):
    """Export hierarchy GeoJSON plus monthly and annual SD JSON payloads."""
    print(f"\n--- Processing hierarchy: {output_prefix} ---")
    if not os.path.exists(shapefile_path):
        print(f"Warning: {shapefile_path} not found. Skipping.")
        return

    try:
        gdf_hierarchy = gpd.read_file(shapefile_path).to_crs(gdf_estaciones_sd.crs)
    except Exception as exc:
        print(f"Critical error reading {shapefile_path}: {exc}")
        return

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
    gdf_filtrado.to_file(geojson_path, driver="GeoJSON")
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

    df_daily = cargar_series_temporales(SD_FILE_PATH)
    if df_daily is None or df_daily.empty:
        print("No valid SD series found. Aborting.")
        return

    station_codes_with_data = set(df_daily.columns)
    gdf_estaciones_sd = gdf_estaciones[gdf_estaciones["code_internal"].isin(station_codes_with_data)].copy()
    gdf_estaciones_sd["variables_disponibles"] = "SD"

    if gdf_estaciones_sd.empty:
        print("No station metadata matched the SD dataset. Aborting.")
        return

    codes_only_in_data = station_codes_with_data - set(gdf_estaciones["code_internal"])
    if codes_only_in_data:
        print(f"Warning: {len(codes_only_in_data)} SD stations are present in the data file but missing in estaciones.csv")
        print(sorted(list(codes_only_in_data))[:10], "..." if len(codes_only_in_data) > 10 else "")

    df_monthly = calcular_agregacion_mensual(df_daily)
    df_annual = calcular_mediana_anual_mayo_octubre(df_daily, df_monthly)

    estaciones_geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    gdf_estaciones_sd[gdf_estaciones_sd.is_valid].to_file(estaciones_geojson_path, driver="GeoJSON")
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
    procesar_jerarquia_geoespacial(PATH_CUENCAS_SHP, ID_COLUMNA_CUENCA, "cuencas", gdf_estaciones_sd, df_monthly, df_annual)
    procesar_jerarquia_geoespacial(PATH_SUBCUENCAS_SHP, ID_COLUMNA_SUBCUENCA, "subcuencas", gdf_estaciones_sd, df_monthly, df_annual)
    procesar_jerarquia_geoespacial(PATH_SUBSUBCUENCAS_SHP, ID_COLUMNA_SUBSUBCUENCA, "subsubcuencas", gdf_estaciones_sd, df_monthly, df_annual)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    version_for_cache = datetime.now().strftime("%Y-%m-%d-sd")
    metadata_path = os.path.join(OUTPUT_DIR, "export_metadata.json")
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "export_timestamp": timestamp,
                "data_version": version_for_cache,
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