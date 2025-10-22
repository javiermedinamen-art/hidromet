import pandas as pd
import geopandas as gpd
import os
import json
import sys
from shapely.geometry import Point

# --------------------------------------------------------------------------
# --- CONFIGURACIÓN DE RUTAS Y COLUMNAS ---
# --------------------------------------------------------------------------
OUTPUT_DIR = "data_static"

# --- CONFIGURACIÓN DE VARIABLES ---
VARIABLES_CONFIG = {
    'ppt': {
        'path': "data/PPT_clean_v1.csv",
        'unit': "mm",
        'display_name': "Precipitacion",
        'fill_rule_percent': 0.8,
        'agg_monthly': 'sum',
        'run_annual': True
    },
    'sd': {
        'path': "data/SD_clean_v2.csv", # Asegúrate que este sea el nombre correcto
        'unit': "cm",
        'display_name': "Altura_de_Nieve",
        'fill_rule_percent': 0.5,
        'agg_monthly': 'mean',
        'run_annual': False
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

# --------------------------------------------------------------------------

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

    # Estandariza todos los code_internal a minúsculas
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
            # Intenta convertir a numérico, reemplazando comas decimales si es necesario
            if df_series[col].dtype == 'object':
                 df_series[col] = df_series[col].astype(str).str.replace(',', '.', regex=False)
            df_series[col] = pd.to_numeric(df_series[col], errors='coerce')

        df_series.dropna(axis=1, how='all', inplace=True)
        print(f"✅ Series cargadas para {variable_name}. Columnas: {len(df_series.columns)}")
        return df_series
    except Exception as e:
        print(f"Error al cargar series temporales de {file_path}: {e}")
        return None

def procesar_jerarquia_geoespacial(shapefile_path, id_columna, output_prefix, gdf_estaciones, df_series, variable_code, var_config):
    """Procesa una capa geoespacial para una variable específica."""
    print(f"\n--- Procesando jerarquía: '{output_prefix}' para variable '{variable_code}' ---")

    fill_rule_percent = var_config['fill_rule_percent']
    agg_method = var_config['agg_monthly']
    run_annual_analysis = var_config['run_annual']

    if not os.path.exists(shapefile_path):
        print(f"Advertencia: No se encontró el shapefile en '{shapefile_path}'. Saltando.")
        return
    try:
        gdf_jerarquia = gpd.read_file(shapefile_path)
        gdf_jerarquia = gdf_jerarquia.to_crs(gdf_estaciones.crs)
    except Exception as e:
        print(f"Error CRÍTICO al leer '{shapefile_path}': {e}")
        return

    codigos_con_variable_actual = df_series.columns
    gdf_estaciones_filtradas = gdf_estaciones[gdf_estaciones['code_internal'].isin(codigos_con_variable_actual)].copy()
    if gdf_estaciones_filtradas.empty:
        print(f"Advertencia: Ninguna estación en '{output_prefix}' tiene datos para '{variable_code}'. Saltando.")
        return

    gdf_estaciones_con_jerarquia = gpd.sjoin(gdf_estaciones_filtradas, gdf_jerarquia, how="inner", predicate="within")
    if gdf_estaciones_con_jerarquia.empty:
        print(f"Advertencia: Ninguna estación CON DATOS de '{variable_code}' encontrada dentro de los polígonos de '{output_prefix}'.")
        return

    ids_con_estaciones = gdf_estaciones_con_jerarquia[id_columna].unique()
    gdf_filtrado = gdf_jerarquia[gdf_jerarquia[id_columna].isin(ids_con_estaciones)].copy()

    station_counts = gdf_estaciones_con_jerarquia.groupby(id_columna).size()
    station_counts.name = 'n_estaciones'
    gdf_filtrado = gdf_filtrado.merge(station_counts, left_on=id_columna, right_index=True, how='left')
    gdf_filtrado['n_estaciones'] = gdf_filtrado['n_estaciones'].fillna(0).astype(int)
    print(f"✅ Se calculó el número de estaciones ({variable_code}) para {len(gdf_filtrado)} polígonos.")

    geojson_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{variable_code}.geojson")
    gdf_filtrado.to_file(geojson_path, driver="GeoJSON")
    print(f"✅ GeoJSON de '{output_prefix}' ({variable_code}) exportado a: {geojson_path}")

    count_poligonos = 0
    for poly_id in ids_con_estaciones:
        try: # Añadido bloque try/except por polígono
            estaciones_en_poligono = gdf_estaciones_con_jerarquia[gdf_estaciones_con_jerarquia[id_columna] == poly_id].copy()
            estaciones_en_poligono['elevation'] = pd.to_numeric(estaciones_en_poligono['elevation'], errors='coerce')
            estaciones_ordenadas = estaciones_en_poligono.sort_values(by='elevation', ascending=False, na_position='last')
            codigos_estaciones_ordenados = estaciones_ordenadas['code_internal'].tolist()

            # Asegúrate de que df_series_poligono solo contenga columnas válidas
            columnas_validas = [col for col in codigos_estaciones_ordenados if col in df_series.columns]
            if not columnas_validas:
                 print(f"  Polígono {poly_id}: No hay columnas válidas en df_series. Saltando.")
                 continue
            df_series_poligono = df_series[columnas_validas] # Reindexa solo con válidas

            monthly_agg = df_series_poligono.resample('MS').agg([agg_method, 'count'])

            for col in columnas_validas: # Itera solo sobre columnas válidas
                if (col, 'count') in monthly_agg.columns:
                    days_in_month = monthly_agg.index.days_in_month.to_numpy() # Convertir a numpy array
                    count_values = monthly_agg[(col, 'count')].to_numpy() # Convertir a numpy array
                    # Comparación vectorial
                    condition = count_values < days_in_month * fill_rule_percent
                    # Asignación usando .loc con la condición numpy
                    monthly_agg.loc[condition, (col, agg_method)] = None

            df_agg_monthly = monthly_agg.xs(agg_method, level=1, axis=1).dropna(how='all')

            if df_agg_monthly.empty:
                print(f"  Polígono {poly_id}: No quedaron datos mensuales válidos tras aplicar reglas. Saltando JSON.")
                # Asegúrate de que no exista un archivo JSON vacío o corrupto de una ejecución anterior
                json_path_check = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}.json")
                if os.path.exists(json_path_check): os.remove(json_path_check)
                annual_json_path_check = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual.json")
                if os.path.exists(annual_json_path_check): os.remove(annual_json_path_check)
                continue # Salta a la siguiente iteración del bucle for poly_id

            start_date = df_agg_monthly.index.min()
            end_date = df_agg_monthly.index.max()
            complete_date_index = pd.date_range(start_date, end_date, freq='MS')
            df_complete = df_agg_monthly.reindex(complete_date_index)

            output_json = {}
            for station_code in columnas_validas: # Itera solo sobre columnas válidas
                if station_code not in df_complete.columns:
                    continue
                # Usa .loc para evitar potential KeyError si el índice no es único (aunque no debería pasar aquí)
                station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                 # Asegúrate que station_info no sea una Serie vacía o tenga múltiples filas si hay duplicados
                if isinstance(station_info, pd.DataFrame):
                    station_info = station_info.iloc[0] # Toma la primera fila si hay duplicados

                values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_complete[station_code]]
                elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                output_json[station_code] = {
                    'name': station_info['name'],
                    'elevation': elevation_value,
                    'dates': df_complete.index.strftime('%Y-%m-%d').tolist(),
                    'values': values_with_nulls
                }

            json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}.json")
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(output_json, f, ensure_ascii=False, indent=2) # Añadido indent para debug

            # --- CÁLCULO ANUAL CONDICIONAL ---
            if run_annual_analysis:
                df_sum_monthly_ppt = df_agg_monthly
                hsm_por_estacion = {}
                # ... (resto de lógica anual idéntica)...
                for station_code in df_sum_monthly_ppt.columns:
                     climatologia_mensual = df_sum_monthly_ppt[station_code].groupby(df_sum_monthly_ppt.index.month).mean()
                     if climatologia_mensual.empty or climatologia_mensual.isnull().all(): continue # Skip if no valid monthly data
                     climatologia_ordenada = climatologia_mensual.sort_values(ascending=False)
                     umbral_precipitacion = climatologia_mensual.sum() * 0.75
                     suma_acumulada = climatologia_ordenada.cumsum()
                     meses_significativos = suma_acumulada[suma_acumulada <= umbral_precipitacion].index.tolist()
                     if not meses_significativos and not climatologia_ordenada.empty:
                         meses_significativos = [climatologia_ordenada.index[0]]
                     hsm_por_estacion[station_code] = meses_significativos

                df_annual = df_sum_monthly_ppt.resample('YS').sum(min_count=1)
                monthly_counts = df_sum_monthly_ppt.resample('YS').count()

                for station_code in df_annual.columns:
                    for year_start_date in df_annual.index:
                        year = year_start_date.year
                        try: num_meses_validos = monthly_counts.loc[year_start_date, station_code]
                        except KeyError: num_meses_validos = 0
                        if num_meses_validos < 9:
                            df_annual.loc[year_start_date, station_code] = None
                            continue
                        meses_clave = hsm_por_estacion.get(station_code, [])
                        if not meses_clave: continue

                        # Asegúrate que los datos del año existan antes de acceder
                        if str(year) not in df_sum_monthly_ppt.index.year.astype(str):
                             df_annual.loc[year_start_date, station_code] = None
                             continue
                        datos_del_ano = df_sum_monthly_ppt.loc[str(year), station_code]
                        if isinstance(datos_del_ano, pd.Series): # Asegura que sea una Serie (no un solo valor)
                             for mes_clave in meses_clave:
                                 # Revisa si el mes existe Y si el valor es nulo
                                 month_data = datos_del_ano[datos_del_ano.index.month == mes_clave]
                                 if month_data.empty or pd.isna(month_data.iloc[0]):
                                     df_annual.loc[year_start_date, station_code] = None
                                     break
                        else: # Si datos_del_ano no es una Serie (caso raro), anula el año
                             df_annual.loc[year_start_date, station_code] = None


                if not df_annual.dropna(how='all').empty: # Solo guarda si hay algún dato anual no nulo
                    output_annual_json = {}
                    for station_code in columnas_validas: # Usa columnas válidas
                        if station_code not in df_annual.columns: continue
                        station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                        if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]

                        values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual[station_code]]
                        elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
                        output_annual_json[station_code] = {
                            'name': station_info['name'], 'elevation': elevation_value,
                            'years': df_annual.index.year.tolist(), 'values': values_with_nulls
                        }

                    annual_json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual.json")
                    with open(annual_json_path, 'w', encoding='utf-8') as f:
                        json.dump(output_annual_json, f, ensure_ascii=False, indent=2) # Añadido indent
                else:
                    print(f"  Polígono {poly_id}: No quedaron datos anuales válidos. Saltando JSON anual.")
                     # Asegúrate de eliminar archivo anual si existe de antes
                    annual_json_path_check = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual.json")
                    if os.path.exists(annual_json_path_check): os.remove(annual_json_path_check)
            else:
                 pass # No se ejecuta análisis anual para esta variable
                 # Asegúrate de eliminar archivo anual si existe de antes
                 annual_json_path_check = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_{variable_code}_annual.json")
                 if os.path.exists(annual_json_path_check): os.remove(annual_json_path_check)

            count_poligonos += 1

        except Exception as e_poly: # Captura error por polígono
             print(f"ERROR procesando polígono {poly_id} para {variable_code}: {e_poly}")
             # Considera continuar con el siguiente polígono
             continue


    print(f"✅ Exportados {count_poligonos} archivos JSON de datos para '{output_prefix}' ({variable_code}).")

def exportar_datos_estaticos():
    """Función principal."""
    print("Iniciando exportación de datos estáticos...")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # PASO 1: Cargar estaciones
    gdf_estaciones = cargar_estaciones()
    if gdf_estaciones is None or gdf_estaciones.empty:
        print("Fallo la carga de estaciones. Abortando.")
        return

    # PASO 2: Generación dinámica de variables disponibles
    print("\n--- Analizando disponibilidad de variables ---")
    station_vars = {}
    all_station_codes_in_data = set() # Para verificar consistencia

    for var_code, config in VARIABLES_CONFIG.items():
        file_path = config['path']
        if not os.path.exists(file_path):
            print(f"Advertencia: No se encontró {file_path}. Saltando {var_code}.")
            continue
        try:
            df_check = pd.read_csv(file_path, nrows=0, encoding='utf-8', sep=',', decimal='.')
            df_check.columns = df_check.columns.astype(str).str.strip().str.lower()
            station_codes_with_data = [col for col in df_check.columns if col != 'date']
            all_station_codes_in_data.update(station_codes_with_data) # Agrega códigos al set global

            for code in station_codes_with_data:
                if code not in station_vars: station_vars[code] = []
                station_vars[code].append(var_code.upper())
            print(f"Detectadas {len(station_codes_with_data)} estaciones para '{var_code}'.")
        except Exception as e:
            print(f"Error leyendo header de {file_path}: {e}")

    # Añade columna 'variables_disponibles'
    def get_vars_string(row):
        vars_list = station_vars.get(row['code_internal'], [])
        return '|'.join(sorted(vars_list)) # Ordena para consistencia

    gdf_estaciones['variables_disponibles'] = gdf_estaciones.apply(get_vars_string, axis=1)
    print("✅ Información de variables añadida a gdf_estaciones.")

    # PASO 3: Exportar GeoJSON y Nombres de estaciones
    print("\n--- Exportando GeoJSON y Nombres ---")
    geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    gdf_estaciones.to_file(geojson_path, driver="GeoJSON")
    print(f"✅ GeoJSON de estaciones exportado a: {geojson_path}")

    try: # Exporta solo nombres de estaciones que están en gdf_estaciones
        station_names = gdf_estaciones.set_index('code_internal')['name'].to_dict()
        with open(os.path.join(OUTPUT_DIR, "station_names.json"), 'w', encoding='utf-8') as f:
            json.dump(station_names, f, ensure_ascii=False, indent=2)
        print("✅ Archivo de nombres de estaciones exportado.")
    except Exception as e:
        print(f"Advertencia: no se pudo exportar station_names.json: {e}")

     # --- Verificación de Consistencia ---
    codes_in_gdf = set(gdf_estaciones['code_internal'])
    codes_only_in_data = all_station_codes_in_data - codes_in_gdf
    if codes_only_in_data:
         print(f"\nADVERTENCIA: Las siguientes estaciones tienen datos CSV pero NO están en '{PATH_ESTACIONES}':")
         print(sorted(list(codes_only_in_data)))
         print("Estos datos no se incluirán en los análisis geoespaciales.")
     # --- Fin Verificación ---

    # PASO 4: Procesar cada variable
    for var_code, config in VARIABLES_CONFIG.items():
        print(f"\n=======================================================")
        print(f"--- Iniciando procesamiento para variable: {var_code.upper()} ---")
        print(f"=======================================================")

        df_series = cargar_series_temporales(config['path'], var_code)
        if df_series is None:
            print(f"Fallo la carga de series para {var_code}. Saltando.")
            continue

        # Exportar CSVs individuales
        print(f"\n--- Exportando CSV individuales para '{var_code}' ---")
        count = 0
        for code in df_series.columns:
            if code not in codes_in_gdf: # Usa el set precalculado
                # Ya advertimos antes, así que solo saltamos
                continue

            series_to_export = df_series[[code]].copy().dropna()
            if series_to_export.empty: continue # Salta si no hay datos no nulos

            col_name = f"{config['display_name']}_{config['unit']}"
            series_to_export.columns = [col_name]
            csv_path = os.path.join(OUTPUT_DIR, f"{code}_{var_code}.csv")
            series_to_export.to_csv(csv_path, index=True, index_label="Fecha")
            count += 1
        print(f"✅ Exportados {count} archivos CSV para '{var_code}'.")

        # Procesar jerarquías
        procesar_jerarquia_geoespacial(PATH_CUENCAS_SHP, ID_COLUMNA_CUENCA, 'cuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBCUENCAS_SHP, ID_COLUMNA_SUBCUENCA, 'subcuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBSUBCUENCAS_SHP, ID_COLUMNA_SUBSUBCUENCA, 'subsubcuencas', gdf_estaciones, df_series, var_code, config)

    print("\n--- Proceso de exportación finalizado ---")


if __name__ == "__main__":
    exportar_datos_estaticos()