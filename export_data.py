import pandas as pd
import geopandas as gpd
import os
import json
import sys
from shapely.geometry import Point
import numpy as np # Necesario para isnan

# --- (CONFIGURACIÓN DE RUTAS Y VARIABLES_CONFIG sin cambios) ---
OUTPUT_DIR = "data_static"
VARIABLES_CONFIG = {
    'ppt': {
        'path': "data/PPT_clean_v1.csv", 'unit': "mm", 'display_name': "Precipitacion",
        'fill_rule_percent': 0.8, 'agg_monthly': 'sum', 'run_annual': True
    },
    'sd': {
        'path': "data/SD_clean_v2.csv", 'unit': "cm", 'display_name': "Altura_de_Nieve",
        'fill_rule_percent': 0.5, 'agg_monthly': 'mean', 'run_annual': False
    },
    'at': {
        'path': "data/AT_clean_v1.csv", 'unit': "°C", 'display_name': "Temperatura",
        'fill_rule_percent': 0.8, 'agg_monthly': 'mean', 'run_annual': False
    }
}
PATH_ESTACIONES = "data/estaciones.csv"
PATH_CUENCAS_SHP = "data/cuencas_bna/cuencas_bna.shp"
ID_COLUMNA_CUENCA = 'COD_CUEN'
NOMBRE_COLUMNA_CUENCA = 'NOM_CUEN'
PATH_SUBCUENCAS_SHP = "data/subcuencas_bna/subcuencas_bna.shp"
ID_COLUMNA_SUBCUENCA = 'COD_SUBC'
NOMBRE_COLUMNA_SUBCUENCA = 'NOM_SUBC'
PATH_SUBSUBCUENCAS_SHP = "data/subsubcuencas_bna/subsubcuencas_bna.shp"
ID_COLUMNA_SUBSUBCUENCA = 'COD_SSUBC'
NOMBRE_COLUMNA_SUBSUBCUENCA = 'NOM_SSUBC'
# --- (FIN CONFIGURACIÓN) ---


# --- (cargar_estaciones y cargar_series_temporales sin cambios) ---
def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    if not os.path.exists(PATH_ESTACIONES):
        print(f"Error CRÍTICO: No se encuentra el archivo de estaciones en {PATH_ESTACIONES}")
        return None
    df = None
    try:
        df = pd.read_csv(PATH_ESTACIONES, encoding="utf-8", sep=',', decimal='.')
        print(f"Leyendo estaciones desde: {PATH_ESTACIONES}")
    except Exception:
        try:
            df = pd.read_csv(PATH_ESTACIONES, encoding="cp1252", sep=',', decimal='.')
            print(f"Leyendo estaciones desde: {PATH_ESTACIONES} con cp1252")
        except Exception as e:
            print(f"Error CRÍTICO al leer {PATH_ESTACIONES}. Error: {e}")
            return None

    df.columns = [c.strip().lower() for c in df.columns]
    required_cols = ['lon', 'lat', 'code_internal', 'name', 'fuente', 'basin']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Faltan columnas requeridas en estaciones.csv. Necesarias: {required_cols}")
        print(f"Columnas encontradas: {list(df.columns)}")
        return None

    df['longitude'] = pd.to_numeric(df['lon'].astype(str).str.strip(), errors='coerce')
    df['latitude'] = pd.to_numeric(df['lat'].astype(str).str.strip(), errors='coerce')
    df_clean = df.dropna(subset=['longitude', 'latitude']).copy()
    df_clean['code_internal'] = df_clean['code_internal'].astype(str).str.strip().str.lower()
    if df_clean.empty:
        print("Advertencia: No quedan estaciones con coordenadas válidas.")
        return gpd.GeoDataFrame()

    geometry = [Point(xy) for xy in zip(df_clean['longitude'], df_clean['latitude'])]
    gdf = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    return gdf

def cargar_series_temporales(file_path, variable_name):
    """Carga los datos de series temporales, asegura la conversión a datetime y numérico."""
    print(f"\n--- Cargando series para variable: '{variable_name}' ---")
    if not os.path.exists(file_path):
        print(f"Error: no se encontró el archivo de series: {file_path}")
        return None
    try:
        df = pd.read_csv(file_path, encoding="utf-8", sep=',', decimal='.')
        df.columns = df.columns.astype(str).str.strip()
        df.rename(columns={c: c.lower() for c in df.columns}, inplace=True)
        if 'date' not in df.columns:
            print(f"Error CRÍTICO: La columna 'date' no se encontró en {file_path}.")
            return None
        df['date'] = pd.to_datetime(df['date'].astype(str).str.strip(), errors='coerce')
        df.dropna(subset=['date'], inplace=True)
        df_series = df.set_index('date')
        for col in df_series.columns:
            if col == 'date': continue
            if df_series[col].dtype == 'object':
                 df_series[col] = df_series[col].astype(str).str.replace(',', '.', regex=False)
            df_series[col] = pd.to_numeric(df_series[col], errors='coerce')
        df_series.dropna(axis=1, how='all', inplace=True)
        df_series.dropna(axis=0, how='all', inplace=True)
        print(f"✅ Series cargadas para {variable_name}. Columnas: {len(df_series.columns)}, Filas: {len(df_series)}")
        if df_series.empty:
             print(f"ADVERTENCIA: El DataFrame de series para {variable_name} está vacío después de cargar y limpiar.")
        return df_series
    except Exception as e:
        print(f"Error al cargar series temporales de {file_path}: {e}")
        return None
# --- (FIN cargar_estaciones y cargar_series_temporales) ---


def procesar_jerarquia_geoespacial(shapefile_path, id_columna, output_prefix, gdf_estaciones, df_series, variable_code, var_config):
    """Procesa una capa geoespacial para una variable específica."""
    print(f"\n--- Procesando jerarquía: '{output_prefix}' para variable '{variable_code}' ---")

    fill_rule_percent = var_config['fill_rule_percent']
    agg_method = var_config['agg_monthly']
    run_annual_analysis_ppt = var_config['run_annual'] # Renombrado para claridad

    # --- (Lógica inicial de carga de shapefile y filtrado de estaciones sin cambios) ---
    if not os.path.exists(shapefile_path): print(f"Advertencia: No se encontró {shapefile_path}. Saltando."); return
    try:
        gdf_jerarquia = gpd.read_file(shapefile_path).to_crs(gdf_estaciones.crs)
    except Exception as e: print(f"Error CRÍTICO al leer '{shapefile_path}': {e}"); return

    codigos_con_variable_actual = df_series.columns
    gdf_estaciones_filtradas = gdf_estaciones[gdf_estaciones['code_internal'].isin(codigos_con_variable_actual)].copy()
    if gdf_estaciones_filtradas.empty: print(f"Advertencia: Ninguna estación en {output_prefix} tiene datos para {variable_code}. Saltando."); return

    gdf_estaciones_con_jerarquia = gpd.sjoin(gdf_estaciones_filtradas, gdf_jerarquia, how="inner", predicate="within")
    if gdf_estaciones_con_jerarquia.empty: print(f"Advertencia: Ninguna estación CON DATOS de {variable_code} encontrada en {output_prefix}."); return

    ids_con_estaciones = gdf_estaciones_con_jerarquia[id_columna].unique()
    gdf_filtrado = gdf_jerarquia[gdf_jerarquia[id_columna].isin(ids_con_estaciones)].copy()

    station_counts = gdf_estaciones_con_jerarquia.groupby(id_columna).size().rename('n_estaciones')
    gdf_filtrado = gdf_filtrado.merge(station_counts, left_on=id_columna, right_index=True, how='left')
    gdf_filtrado['n_estaciones'] = gdf_filtrado['n_estaciones'].fillna(0).astype(int)
    print(f"✅ {len(gdf_filtrado)} polígonos con estaciones ({variable_code}).")

    geojson_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{variable_code}.geojson")
    try:
        gdf_filtrado.to_file(geojson_path, driver="GeoJSON")
        print(f"✅ GeoJSON exportado: {geojson_path}")
    except Exception as e_geojson: print(f"ERROR guardando GeoJSON {geojson_path}: {e_geojson}"); return
    # --- (FIN Lógica inicial) ---

    count_poligonos = 0
    for poly_id in ids_con_estaciones:
        # Define paths para eliminar en caso de error o datos vacíos
        json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}.json")
        annual_ppt_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual.json") # Solo PPT
        annual_sd_max_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual_max.json") # Solo SD
        annual_at_stats_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual_stats.json") # Solo AT

        try:
            estaciones_en_poligono = gdf_estaciones_con_jerarquia[gdf_estaciones_con_jerarquia[id_columna] == poly_id].copy()
            estaciones_en_poligono['elevation'] = pd.to_numeric(estaciones_en_poligono['elevation'], errors='coerce')
            estaciones_ordenadas = estaciones_en_poligono.sort_values(by='elevation', ascending=False, na_position='last')
            codigos_estaciones_ordenados = estaciones_ordenadas['code_internal'].tolist()

            columnas_validas = [col for col in codigos_estaciones_ordenados if col in df_series.columns]
            if not columnas_validas: print(f"  Polígono {poly_id}: Sin columnas válidas. Saltando."); continue
            df_series_poligono = df_series[columnas_validas]

            monthly_agg = df_series_poligono.resample('MS').agg([agg_method, 'count'])

            for col in columnas_validas:
                if (col, 'count') in monthly_agg.columns:
                    days_in_month = monthly_agg.index.map(lambda dt: dt.daysinmonth)
                    count_values = monthly_agg[(col, 'count')].values
                    # Asegurar que days_in_month tenga la misma longitud que count_values
                    if len(days_in_month) != len(count_values):
                         # Esto puede pasar si hay NaNs al principio/final que se pierden en .values
                         # Reindexar count_values para alinear con el índice de monthly_agg
                         count_series = monthly_agg[(col, 'count')]
                         condition = count_series < days_in_month.to_series(index=count_series.index) * fill_rule_percent
                    else:
                         condition = count_values < days_in_month * fill_rule_percent

                    monthly_agg.loc[condition, (col, agg_method)] = None

            df_agg_monthly = monthly_agg.xs(agg_method, level=1, axis=1).dropna(how='all')

            if df_agg_monthly.empty:
                print(f"  Polígono {poly_id}: Sin datos mensuales válidos. Saltando JSON.")
                # Elimina archivos existentes si no hay datos válidos
                if os.path.exists(json_path): os.remove(json_path)
                if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)
                if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)
                continue

            # --- Guardar JSON Mensual (sin cambios) ---
            start_date = df_agg_monthly.index.min()
            end_date = df_agg_monthly.index.max()
            complete_date_index = pd.date_range(start_date, end_date, freq='MS')
            df_complete = df_agg_monthly.reindex(complete_date_index)
            output_json = {}
            for station_code in columnas_validas:
                if station_code not in df_complete.columns: continue
                try:
                     station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                     if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                except KeyError: station_info = pd.Series({'name': station_code, 'elevation': None})
                values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_complete[station_code]]
                elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                output_json[station_code] = {
                    'name': station_info['name'], 'elevation': elevation_value,
                    'dates': df_complete.index.strftime('%Y-%m-%d').tolist(), 'values': values_with_nulls
                }
            with open(json_path, 'w', encoding='utf-8') as f: json.dump(output_json, f, ensure_ascii=False, indent=2)
            # --- FIN Guardar JSON Mensual ---

            # --- ANÁLISIS ANUAL (CONDICIONAL) ---
            monthly_counts_annual = df_agg_monthly.resample('YS').count() # Conteo de meses válidos por año

            # A) Precipitación (PPT): Análisis robusto original
            if run_annual_analysis_ppt and variable_code == 'ppt':
                print(f"  Polígono {poly_id}: Ejecutando análisis anual robusto PPT...")
                # ... (Lógica idéntica al código anterior para PPT: hsm, df_annual, etc.) ...
                df_sum_monthly_ppt = df_agg_monthly
                hsm_por_estacion = {}
                for station_code in df_sum_monthly_ppt.columns:
                    climatologia_mensual = df_sum_monthly_ppt[station_code].groupby(df_sum_monthly_ppt.index.month).mean()
                    if climatologia_mensual.empty or climatologia_mensual.isnull().all(): continue
                    climatologia_ordenada = climatologia_mensual.sort_values(ascending=False)
                    umbral_precipitacion = climatologia_mensual.sum() * 0.75
                    suma_acumulada = climatologia_ordenada.cumsum()
                    meses_significativos = suma_acumulada[suma_acumulada <= umbral_precipitacion].index.tolist()
                    if not meses_significativos and not climatologia_ordenada.empty: meses_significativos = [climatologia_ordenada.index[0]]
                    hsm_por_estacion[station_code] = meses_significativos

                df_annual = df_sum_monthly_ppt.resample('YS').sum(min_count=1)

                for station_code in df_annual.columns:
                    for year_start_date in df_annual.index:
                        year = year_start_date.year
                        try: num_meses_validos = monthly_counts_annual.loc[year_start_date, station_code]
                        except KeyError: num_meses_validos = 0
                        if num_meses_validos < 9: df_annual.loc[year_start_date, station_code] = None; continue
                        meses_clave = hsm_por_estacion.get(station_code, [])
                        if not meses_clave: continue
                        if str(year) not in df_sum_monthly_ppt.index.year.astype(str): df_annual.loc[year_start_date, station_code] = None; continue
                        datos_del_ano_series = df_sum_monthly_ppt.loc[str(year), station_code]
                        if not isinstance(datos_del_ano_series, pd.Series):
                            if pd.api.types.is_scalar(datos_del_ano_series):
                                original_index = df_sum_monthly_ppt.loc[str(year)].index
                                if not original_index.empty: datos_del_ano_series = pd.Series([datos_del_ano_series], index=[original_index[0]])
                                else: df_annual.loc[year_start_date, station_code] = None; continue
                            else: df_annual.loc[year_start_date, station_code] = None; continue
                        for mes_clave in meses_clave:
                             month_data = datos_del_ano_series[datos_del_ano_series.index.month == mes_clave]
                             if month_data.empty or pd.isna(month_data.iloc[0]): df_annual.loc[year_start_date, station_code] = None; break

                if not df_annual.dropna(how='all').empty:
                    output_annual_json = {}
                    # ... (creación de output_annual_json idéntica) ...
                    for station_code in columnas_validas:
                         if station_code not in df_annual.columns: continue
                         try:
                             station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                             if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                         except KeyError: station_info = pd.Series({'name': station_code, 'elevation': None})
                         values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual[station_code]]
                         elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                         output_annual_json[station_code] = {'name': station_info['name'], 'elevation': elevation_value, 'years': df_annual.index.year.tolist(), 'values': values_with_nulls}

                    with open(annual_ppt_path, 'w', encoding='utf-8') as f: json.dump(output_annual_json, f, ensure_ascii=False, indent=2)
                    print(f"    ✅ Guardado: {annual_ppt_path}")
                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales PPT válidos.")
                    if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                # Limpia los otros archivos anuales si existen por error
                if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)
                if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)


            # B) Altura de Nieve (SD): Máximo mensual anual
            elif variable_code == 'sd':
                print(f"  Polígono {poly_id}: Calculando máximo anual SD...")
                df_annual_max = df_agg_monthly.resample('YS').max() # Calcula el máximo mensual por año

                # Aplica regla de validez: Nulo si menos de 6 meses válidos
                for station_code in df_annual_max.columns:
                     for year_start_date in df_annual_max.index:
                         try: num_meses_validos = monthly_counts_annual.loc[year_start_date, station_code]
                         except KeyError: num_meses_validos = 0
                         if num_meses_validos < 6: # ¡Regla de 6 meses para SD max!
                             df_annual_max.loc[year_start_date, station_code] = None

                if not df_annual_max.dropna(how='all').empty:
                    output_annual_max_json = {}
                    for station_code in columnas_validas:
                        if station_code not in df_annual_max.columns: continue
                        try:
                            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                            if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                        except KeyError: station_info = pd.Series({'name': station_code, 'elevation': None})

                        # Asegura que los valores NaN de numpy se conviertan a None
                        values_with_nulls = [None if np.isnan(val) else round(val, 1) for val in df_annual_max[station_code].fillna(np.nan)]
                        elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                        output_annual_max_json[station_code] = {
                            'name': station_info['name'], 'elevation': elevation_value,
                            'years': df_annual_max.index.year.tolist(), 'values': values_with_nulls # 'values' contiene los máximos
                        }

                    with open(annual_sd_max_path, 'w', encoding='utf-8') as f: json.dump(output_annual_max_json, f, ensure_ascii=False, indent=2)
                    print(f"    ✅ Guardado: {annual_sd_max_path}")
                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales SD max válidos.")
                    if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)
                # Limpia los otros archivos anuales
                if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)


            # C) Temperatura (AT): Máximo y Mínimo mensual anual
            elif variable_code == 'at':
                print(f"  Polígono {poly_id}: Calculando extremos anuales AT...")
                df_annual_stats = df_agg_monthly.resample('YS').agg(['max', 'min']) # Calcula max y min

                # Aplica regla de validez: Nulo si menos de 9 meses válidos
                for station_code in df_annual_stats.columns.levels[0]: # Itera sobre códigos de estación
                     for year_start_date in df_annual_stats.index:
                         try: num_meses_validos = monthly_counts_annual.loc[year_start_date, station_code]
                         except KeyError: num_meses_validos = 0
                         if num_meses_validos < 9: # ¡Regla de 9 meses para AT stats!
                             df_annual_stats.loc[year_start_date, (station_code, 'max')] = None
                             df_annual_stats.loc[year_start_date, (station_code, 'min')] = None

                # Verifica si hay *algún* dato válido (max o min) antes de guardar
                if not df_annual_stats.dropna(how='all').empty:
                    output_annual_stats_json = {}
                    for station_code in columnas_validas:
                        if station_code not in df_annual_stats.columns.levels[0]: continue
                        try:
                            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                            if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                        except KeyError: station_info = pd.Series({'name': station_code, 'elevation': None})

                        # Extrae max y min, convirtiendo NaN a None
                        max_values = [None if np.isnan(val) else round(val, 1) for val in df_annual_stats[(station_code, 'max')].fillna(np.nan)]
                        min_values = [None if np.isnan(val) else round(val, 1) for val in df_annual_stats[(station_code, 'min')].fillna(np.nan)]

                        # Solo incluye la estación si tiene al menos un max o min no nulo
                        if any(v is not None for v in max_values) or any(v is not None for v in min_values):
                             elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                             output_annual_stats_json[station_code] = {
                                'name': station_info['name'], 'elevation': elevation_value,
                                'years': df_annual_stats.index.year.tolist(),
                                'max_values': max_values, # Lista de máximos
                                'min_values': min_values  # Lista de mínimos
                             }

                    # Guarda el archivo solo si el diccionario no está vacío
                    if output_annual_stats_json:
                        with open(annual_at_stats_path, 'w', encoding='utf-8') as f: json.dump(output_annual_stats_json, f, ensure_ascii=False, indent=2)
                        print(f"    ✅ Guardado: {annual_at_stats_path}")
                    else:
                        print(f"  Polígono {poly_id}: Sin datos anuales AT stats válidos (max/min).")
                        if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)

                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales AT stats válidos (max/min).")
                    if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)
                # Limpia los otros archivos anuales
                if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)

            else: # Caso para otras variables futuras sin análisis anual
                 if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                 if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)
                 if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)

            count_poligonos += 1

        except Exception as e_poly:
             print(f"ERROR procesando polígono {poly_id} para {variable_code}: {e_poly}")
             # Elimina archivos potencialmente corruptos
             if os.path.exists(json_path): os.remove(json_path)
             if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
             if os.path.exists(annual_sd_max_path): os.remove(annual_sd_max_path)
             if os.path.exists(annual_at_stats_path): os.remove(annual_at_stats_path)
             continue

    print(f"✅ Exportados {count_poligonos} archivos JSON mensuales para '{output_prefix}' ({variable_code}).")


def exportar_datos_estaticos():
    """Función principal."""
    print("Iniciando exportación de datos estáticos...")
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

    gdf_estaciones = cargar_estaciones()
    if gdf_estaciones is None or gdf_estaciones.empty: print("Fallo la carga de estaciones. Abortando."); return

    print("\n--- Analizando disponibilidad de variables ---")
    station_vars = {}
    all_station_codes_in_data = set()
    for var_code, config in VARIABLES_CONFIG.items():
        file_path = config['path']
        if not os.path.exists(file_path): print(f"Advertencia: No se encontró {file_path}. Saltando {var_code}."); continue
        try:
            df_check = pd.read_csv(file_path, nrows=0, encoding='utf-8', sep=',', decimal='.')
        except Exception:
             try: df_check = pd.read_csv(file_path, nrows=0, encoding='cp1252', sep=',', decimal='.')
             except Exception as e_read: print(f"Error leyendo header de {file_path}: {e_read}"); continue
        df_check.columns = df_check.columns.astype(str).str.strip().str.lower()
        station_codes_with_data = [col for col in df_check.columns if col != 'date']
        all_station_codes_in_data.update(station_codes_with_data)
        for code in station_codes_with_data:
            if code not in station_vars: station_vars[code] = []
            station_vars[code].append(var_code.upper())
        print(f"Detectadas {len(station_codes_with_data)} estaciones para '{var_code}'.")

    def get_vars_string(row): return '|'.join(sorted(station_vars.get(row['code_internal'], [])))
    gdf_estaciones['variables_disponibles'] = gdf_estaciones.apply(get_vars_string, axis=1)
    print("✅ Información de variables añadida a gdf_estaciones.")

    print("\n--- Exportando GeoJSON y Nombres ---")
    geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    try:
        gdf_estaciones.to_file(geojson_path, driver="GeoJSON")
        print(f"✅ GeoJSON de estaciones exportado: {geojson_path}")
    except Exception as e_geojson_main: print(f"ERROR CRÍTICO al guardar {geojson_path}: {e_geojson_main}"); return

    try:
        gdf_unique = gdf_estaciones.drop_duplicates(subset=['code_internal'], keep='first')
        # Filtra el gdf_unique ANTES de convertir a dict para asegurar que los nombres existan
        station_names = gdf_unique[gdf_unique['code_internal'].isin(all_station_codes_in_data)].set_index('code_internal')['name'].to_dict()
        with open(os.path.join(OUTPUT_DIR, "station_names.json"), 'w', encoding='utf-8') as f: json.dump(station_names, f, ensure_ascii=False, indent=2)
        print("✅ Archivo de nombres de estaciones exportado.")
    except Exception as e: print(f"Advertencia: no se pudo exportar station_names.json: {e}")

    codes_in_gdf = set(gdf_estaciones['code_internal'])
    codes_only_in_data = all_station_codes_in_data - codes_in_gdf
    if codes_only_in_data:
         print(f"\nADVERTENCIA: {len(codes_only_in_data)} estaciones con datos CSV NO están en '{PATH_ESTACIONES}':")
         print(sorted(list(codes_only_in_data))[:10], "..." if len(codes_only_in_data) > 10 else "")
         print("Estos datos NO se incluirán en análisis geoespaciales ni aparecerán en el mapa.")

    for var_code, config in VARIABLES_CONFIG.items():
        print(f"\n=======================================================")
        print(f"--- Iniciando procesamiento para variable: {var_code.upper()} ---")
        print(f"=======================================================")

        df_series = cargar_series_temporales(config['path'], var_code)
        if df_series is None or df_series.empty: print(f"No hay datos de series válidos para {var_code}. Saltando."); continue

        print(f"\n--- Exportando CSV individuales para '{var_code}' ---")
        count = 0
        stations_processed = set()
        for code in df_series.columns:
            if code in stations_processed or code not in codes_in_gdf: continue
            series_to_export = df_series[[code]].copy().dropna()
            if series_to_export.empty: continue
            col_name = f"{config['display_name']}_{config['unit']}"
            series_to_export.columns = [col_name]
            csv_path = os.path.join(OUTPUT_DIR, f"{code}_{var_code}.csv")
            try:
                 series_to_export.to_csv(csv_path, index=True, index_label="Fecha")
                 count += 1
                 stations_processed.add(code)
            except Exception as e_csv: print(f"ERROR al guardar CSV {csv_path}: {e_csv}")
        print(f"✅ Exportados {count} archivos CSV para '{var_code}'.")

        procesar_jerarquia_geoespacial(PATH_CUENCAS_SHP, ID_COLUMNA_CUENCA, 'cuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBCUENCAS_SHP, ID_COLUMNA_SUBCUENCA, 'subcuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBSUBCUENCAS_SHP, ID_COLUMNA_SUBSUBCUENCA, 'subsubcuencas', gdf_estaciones, df_series, var_code, config)

    print("\n--- Proceso de exportación finalizado ---")

if __name__ == "__main__":
    exportar_datos_estaticos()