import pandas as pd
import geopandas as gpd
import os
import json
import sys
from shapely.geometry import Point
import numpy as np # Necesario para isnan

# --- (CONFIGURACIÓN DE RUTAS Y VARIABLES_CONFIG) ---
OUTPUT_DIR = "data_static"
VARIABLES_CONFIG = {
    'ppt': {
        'path': "data/PPT_clean_v1.csv", 'unit': "mm", 'display_name': "Precipitacion",
        'fill_rule_percent': 0.7, 'agg_monthly': 'sum', 'run_annual': True
    },
    'sd': {
        'path': "data/SD_clean_v2.csv", 'unit': "cm", 'display_name': "Altura_de_Nieve",
        'fill_rule_percent': 0.7, 'agg_monthly': 'median', 'run_annual': False # run_annual no se usa directamente para sd/at aquí
    },
    'at': {
        'path': "data/AT_clean_v1.csv", 'unit': "°C", 'display_name': "Temperatura",
        'fill_rule_percent': 0.7, 'agg_monthly': 'median', 'run_annual': False # run_annual no se usa directamente para sd/at aquí
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
                # Intentar reemplazar comas y convertir
                try:
                    df_series[col] = df_series[col].astype(str).str.replace(',', '.', regex=False)
                    df_series[col] = pd.to_numeric(df_series[col], errors='coerce')
                except Exception as e_conv:
                    print(f"Advertencia: No se pudo convertir la columna '{col}' a numérico en {file_path}. Error: {e_conv}. Se dejará como NaN donde falle.")
                    # Forzar conversión, dejando NaNs donde falle
                    df_series[col] = pd.to_numeric(df_series[col].astype(str).str.replace(',', '.', regex=False), errors='coerce')

            elif not pd.api.types.is_numeric_dtype(df_series[col]):
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
    run_annual_analysis_ppt = var_config.get('run_annual', False) # Usar .get() para compatibilidad

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
        annual_ppt_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_ppt_annual.json") # Solo PPT, nombre específico
        annual_sd_median_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_sd_annual_median.json") # Nombre específico SD
        annual_at_median_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_at_annual_median.json") # Nombre específico AT


        try:
            estaciones_en_poligono = gdf_estaciones_con_jerarquia[gdf_estaciones_con_jerarquia[id_columna] == poly_id].copy()
            estaciones_en_poligono['elevation'] = pd.to_numeric(estaciones_en_poligono['elevation'], errors='coerce')
            estaciones_ordenadas = estaciones_en_poligono.sort_values(by='elevation', ascending=False, na_position='last')
            codigos_estaciones_ordenados = estaciones_ordenadas['code_internal'].tolist()

            columnas_validas = [col for col in codigos_estaciones_ordenados if col in df_series.columns]
            if not columnas_validas: print(f"  Polígono {poly_id}: Sin columnas válidas. Saltando."); continue
            
            df_series_poligono = df_series[columnas_validas] 

            # Asegurar que el índice sea DatetimeIndex
            if not isinstance(df_series_poligono.index, pd.DatetimeIndex):
                df_series_poligono.index = pd.to_datetime(df_series_poligono.index)

            # Recalcular agregación mensual aquí dentro del bucle
            monthly_agg = df_series_poligono.resample('MS').agg([agg_method, 'count'])

            # Aplicar regla de llenado
            for col in columnas_validas:
                if (col, 'count') in monthly_agg.columns:
                    # Usar map para obtener daysinmonth de forma segura
                    days_in_month_series = pd.Series(monthly_agg.index.map(lambda dt: dt.daysinmonth), index=monthly_agg.index)
                    count_series = monthly_agg[(col, 'count')]
                    
                    # Comparación robusta usando índices alineados
                    condition = count_series < days_in_month_series * fill_rule_percent
                    
                    # Aplicar la condición para anular el valor agregado si no se cumple
                    monthly_agg.loc[condition, (col, agg_method)] = np.nan # Usar np.nan es más seguro


            df_agg_monthly = monthly_agg.xs(agg_method, level=1, axis=1).dropna(how='all', axis=0) # Drop rows where all stations are NaN

            if df_agg_monthly.empty:
                print(f"  Polígono {poly_id}: Sin datos mensuales válidos después de aplicar reglas. Saltando JSON.")
                # Elimina archivos existentes si no hay datos válidos
                if os.path.exists(json_path): os.remove(json_path)
                if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                if os.path.exists(annual_sd_median_path): os.remove(annual_sd_median_path) 
                if os.path.exists(annual_at_median_path): os.remove(annual_at_median_path) 
                continue

            # --- Guardar JSON Mensual ---
            start_date = df_agg_monthly.index.min()
            end_date = df_agg_monthly.index.max()
            complete_date_index = pd.date_range(start_date, end_date, freq='MS')
            # Reindexar ANTES de iterar para asegurar índice completo
            df_complete = df_agg_monthly.reindex(complete_date_index) 
            output_json = {}
            for station_code in columnas_validas:
                # Verificar si la estación existe en el df_complete reindexado
                if station_code not in df_complete.columns: continue 
                
                try:
                    # Usar .get() para evitar KeyError si la elevación es NaN o falta
                    station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                    if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                    elevation_value = station_info.get('elevation') # Usar .get()
                    if pd.isna(elevation_value): elevation_value = None # Convertir NaN a None explícitamente
                    station_name = station_info.get('name', station_code) # Usar código si falta nombre
                except KeyError: 
                    station_name = station_code
                    elevation_value = None
                
                # Asegurar que los valores NaN de numpy/pandas se conviertan a None para JSON
                values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_complete[station_code]]
                
                # Solo añadir la estación si tiene al menos un valor no nulo
                if any(v is not None for v in values_with_nulls):
                    output_json[station_code] = {
                        'name': station_name, 
                        'elevation': elevation_value,
                        'dates': df_complete.index.strftime('%Y-%m-%d').tolist(), 
                        'values': values_with_nulls
                    }
            
            # Guardar solo si el JSON tiene contenido
            if output_json:
                with open(json_path, 'w', encoding='utf-8') as f: json.dump(output_json, f, ensure_ascii=False, indent=2)
            else:
                print(f"  Polígono {poly_id}: Sin datos mensuales válidos para ninguna estación. No se guarda JSON mensual.")
                if os.path.exists(json_path): os.remove(json_path) # Limpiar si existe

            # --- FIN Guardar JSON Mensual ---
            
            # --- ANÁLISIS ANUAL (CONDICIONAL) ---
            # Recalcular conteos anuales basados en df_agg_monthly (que ya tiene la regla aplicada)
            monthly_counts_annual = df_agg_monthly.resample('YS').count() 
            meses_requeridos = [5, 6, 7, 8, 9, 10] # Mayo a Octubre

            # A) Precipitación (PPT): Acumulado anual condicionado a Mayo-Octubre
            if run_annual_analysis_ppt and variable_code == 'ppt':
                print(f"  Polígono {poly_id}: Ejecutando análisis anual PPT (May-Oct)...")
                
                df_annual = df_agg_monthly.resample('YS').sum(min_count=1) 

                for station_code in df_annual.columns:
                    for year_start_date in df_annual.index:
                        year = year_start_date.year
                        
                        try:
                            datos_mensuales_del_ano = df_agg_monthly.loc[str(year), station_code]
                            # Asegurar que sea Series, incluso si solo hay un mes
                            if pd.api.types.is_scalar(datos_mensuales_del_ano):
                                original_index = df_agg_monthly.loc[str(year)].index
                                if not original_index.empty:
                                    datos_mensuales_del_ano = pd.Series([datos_mensuales_del_ano], index=[original_index[0]])
                                else: # Caso raro: año existe pero sin índice?
                                    df_annual.loc[year_start_date, station_code] = np.nan; continue
                            elif not isinstance(datos_mensuales_del_ano, pd.Series): # Si es DataFrame u otro tipo
                                df_annual.loc[year_start_date, station_code] = np.nan; continue
                                
                        except KeyError: 
                            df_annual.loc[year_start_date, station_code] = np.nan; continue

                        # Chequear meses requeridos usando el índice de la serie (que ya tiene NaNs donde no se cumple regla diaria)
                        meses_validos_en_el_ano = datos_mensuales_del_ano.dropna().index.month.unique()
                        condicion_cumplida = all(mes in meses_validos_en_el_ano for mes in meses_requeridos)
                        
                        if not condicion_cumplida:
                            df_annual.loc[year_start_date, station_code] = np.nan 
                
                # Guardar el JSON anual PPT
                if not df_annual.dropna(how='all').empty:
                    output_annual_json = {}
                    for station_code in df_annual.columns: # Iterar sobre columnas del resultado anual
                        # Asegurar que la estación esté en columnas_validas originales
                        if station_code not in columnas_validas: continue 
                        try:
                            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                            if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                            elevation_value = station_info.get('elevation')
                            if pd.isna(elevation_value): elevation_value = None
                            station_name = station_info.get('name', station_code)
                        except KeyError: 
                            station_name = station_code
                            elevation_value = None
                        
                        values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual[station_code]]
                        # Incluir solo si hay algún valor anual válido
                        if any(v is not None for v in values_with_nulls):
                             output_annual_json[station_code] = {'name': station_name, 'elevation': elevation_value, 'years': df_annual.index.year.tolist(), 'values': values_with_nulls}
                    
                    if output_annual_json: # Guardar solo si no está vacío
                        with open(annual_ppt_path, 'w', encoding='utf-8') as f: json.dump(output_annual_json, f, ensure_ascii=False, indent=2)
                        print(f"    ✅ Guardado: {annual_ppt_path}")
                    else:
                        print(f"  Polígono {poly_id}: Sin datos anuales PPT válidos (May-Oct) después del filtro.")
                        if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)

                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales PPT válidos (May-Oct).")
                    if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
                
                # --- CORRECCIÓN: Se eliminaron 2 líneas de limpieza cruzada ---
                # (Líneas que borraban annual_sd_median_path y annual_at_median_path)


            # B) Altura de Nieve (SD): Mediana diaria (May-Oct)
            elif variable_code == 'sd':
                print(f"  Polígono {poly_id}: Calculando mediana diaria (May-Oct) SD...")
                
                df_annual_median = pd.DataFrame(index=monthly_counts_annual.index, columns=df_series_poligono.columns, dtype=float)

                for station_code in df_annual_median.columns:
                    for year_start_date in df_annual_median.index:
                        year = year_start_date.year
                        
                        try:
                            datos_mensuales_del_ano = df_agg_monthly.loc[str(year), station_code]
                            if pd.api.types.is_scalar(datos_mensuales_del_ano):
                                original_index = df_agg_monthly.loc[str(year)].index
                                if not original_index.empty:
                                    datos_mensuales_del_ano = pd.Series([datos_mensuales_del_ano], index=[original_index[0]])
                                else: continue 
                            elif not isinstance(datos_mensuales_del_ano, pd.Series): continue
                        except KeyError: continue 

                        meses_validos_en_el_ano = datos_mensuales_del_ano.dropna().index.month.unique()
                        condicion_mensual_cumplida = all(mes in meses_validos_en_el_ano for mes in meses_requeridos)
                        
                        if not condicion_mensual_cumplida: continue 

                        try:
                            start_period = f"{year}-05-01"
                            end_period = f"{year}-10-31"
                            # Seleccionar datos diarios DENTRO del período May-Oct
                            daily_data_periodo = df_series_poligono.loc[start_period:end_period, station_code]
                            
                            # Calcular mediana SOLO si hay datos no-NaN en el período
                            if not daily_data_periodo.dropna().empty:
                                median_value = daily_data_periodo.median() # median() ignora NaNs pero incluye 0s
                                df_annual_median.loc[year_start_date, station_code] = median_value
                            # Si no hay datos válidos en el período, queda NaN por defecto
                            
                        except Exception as e_median:
                            print(f"    Error calculando mediana SD para {station_code} año {year}: {e_median}")
                
                # Guardar el JSON anual SD
                if not df_annual_median.dropna(how='all').empty:
                    output_annual_median_json = {}
                    for station_code in df_annual_median.columns:
                        if station_code not in columnas_validas: continue
                        try:
                            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                            if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                            elevation_value = station_info.get('elevation')
                            if pd.isna(elevation_value): elevation_value = None
                            station_name = station_info.get('name', station_code)
                        except KeyError: 
                             station_name = station_code
                             elevation_value = None

                        values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual_median[station_code]]
                        if any(v is not None for v in values_with_nulls):
                            output_annual_median_json[station_code] = {
                                'name': station_name, 'elevation': elevation_value,
                                'years': df_annual_median.index.year.tolist(), 
                                'values': values_with_nulls
                            }
                    
                    if output_annual_median_json:
                        # *** ¡ASEGURARSE DE USAR EL PATH CORRECTO AQUÍ! ***
                        with open(annual_sd_median_path, 'w', encoding='utf-8') as f: json.dump(output_annual_median_json, f, ensure_ascii=False, indent=2)
                        print(f"    ✅ Guardado: {annual_sd_median_path}")
                    else:
                         print(f"  Polígono {poly_id}: Sin datos anuales SD median válidos después del filtro.")
                         if os.path.exists(annual_sd_median_path): os.remove(annual_sd_median_path)
                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales SD median válidos.")
                    if os.path.exists(annual_sd_median_path): os.remove(annual_sd_median_path)
                
                # --- CORRECCIÓN: Se eliminaron 2 líneas de limpieza cruzada ---
                # (Líneas que borraban annual_ppt_path y annual_at_median_path)


            # C) Temperatura (AT): Mediana diaria (May-Oct)
            elif variable_code == 'at':
                print(f"  Polígono {poly_id}: Calculando mediana diaria (May-Oct) AT...")
                
                df_annual_median = pd.DataFrame(index=monthly_counts_annual.index, columns=df_series_poligono.columns, dtype=float)

                for station_code in df_annual_median.columns:
                    for year_start_date in df_annual_median.index:
                        year = year_start_date.year
                        
                        try:
                            datos_mensuales_del_ano = df_agg_monthly.loc[str(year), station_code]
                            if pd.api.types.is_scalar(datos_mensuales_del_ano):
                                original_index = df_agg_monthly.loc[str(year)].index
                                if not original_index.empty:
                                    datos_mensuales_del_ano = pd.Series([datos_mensuales_del_ano], index=[original_index[0]])
                                else: continue 
                            elif not isinstance(datos_mensuales_del_ano, pd.Series): continue
                        except KeyError: continue 

                        meses_validos_en_el_ano = datos_mensuales_del_ano.dropna().index.month.unique()
                        condicion_mensual_cumplida = all(mes in meses_validos_en_el_ano for mes in meses_requeridos)
                        
                        if not condicion_mensual_cumplida: continue 

                        try:
                            start_period = f"{year}-05-01"
                            end_period = f"{year}-10-31"
                            daily_data_periodo = df_series_poligono.loc[start_period:end_period, station_code]
                            
                            if not daily_data_periodo.dropna().empty:
                                median_value = daily_data_periodo.median()
                                df_annual_median.loc[year_start_date, station_code] = median_value
                            
                        except Exception as e_median:
                            print(f"    Error calculando mediana AT para {station_code} año {year}: {e_median}")
                            
                # Guardar el JSON anual AT
                if not df_annual_median.dropna(how='all').empty:
                    output_annual_median_json = {}
                    for station_code in df_annual_median.columns:
                         if station_code not in columnas_validas: continue
                         try:
                            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                            if isinstance(station_info, pd.DataFrame): station_info = station_info.iloc[0]
                            elevation_value = station_info.get('elevation')
                            if pd.isna(elevation_value): elevation_value = None
                            station_name = station_info.get('name', station_code)
                         except KeyError: 
                            station_name = station_code
                            elevation_value = None

                         values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual_median[station_code]]
                         if any(v is not None for v in values_with_nulls):
                            output_annual_median_json[station_code] = {
                                'name': station_name, 'elevation': elevation_value,
                                'years': df_annual_median.index.year.tolist(), 
                                'values': values_with_nulls 
                            }
                    
                    if output_annual_median_json:
                        # *** ¡ASEGURARSE DE USAR EL PATH CORRECTO AQUÍ! ***
                        with open(annual_at_median_path, 'w', encoding='utf-8') as f: json.dump(output_annual_median_json, f, ensure_ascii=False, indent=2)
                        print(f"    ✅ Guardado: {annual_at_median_path}")
                    else:
                        print(f"  Polígono {poly_id}: Sin datos anuales AT median válidos después del filtro.")
                        if os.path.exists(annual_at_median_path): os.remove(annual_at_median_path)

                else:
                    print(f"  Polígono {poly_id}: Sin datos anuales AT median válidos.")
                    if os.path.exists(annual_at_median_path): os.remove(annual_at_median_path)
                
                # --- CORRECCIÓN: Se eliminaron 2 líneas de limpieza cruzada ---
                # (Líneas que borraban annual_ppt_path y annual_sd_median_path)
            
            else: # Caso para otras variables futuras sin análisis anual
                # --- CORRECCIÓN: Se eliminaron 3 líneas de limpieza cruzada ---
                # (Líneas que borraban annual_ppt_path, annual_sd_median_path y annual_at_median_path)
                pass # No hacer nada

            # --- FIN DEL BLOQUE ANUAL ---

            count_poligonos += 1

        except Exception as e_poly:
            print(f"ERROR procesando polígono {poly_id} para {variable_code}: {e_poly}")
            # Elimina archivos potencialmente corruptos
            if os.path.exists(json_path): os.remove(json_path)
            if os.path.exists(annual_ppt_path): os.remove(annual_ppt_path)
            if os.path.exists(annual_sd_median_path): os.remove(annual_sd_median_path)
            if os.path.exists(annual_at_median_path): os.remove(annual_at_median_path)
            continue # Continuar con el siguiente polígono

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
    
    # --- INICIO DE CÓDIGO DE "PARTE 2" ---
    
    print("\n--- Exportando GeoJSON y Nombres ---")
    geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    try:
        # Asegurarse de que no haya geometrías inválidas
        gdf_estaciones_valid = gdf_estaciones[gdf_estaciones.is_valid].copy()
        if len(gdf_estaciones_valid) < len(gdf_estaciones):
            print(f"Advertencia: Se excluyeron {len(gdf_estaciones) - len(gdf_estaciones_valid)} estaciones con geometría inválida del GeoJSON.")
        gdf_estaciones_valid.to_file(geojson_path, driver="GeoJSON")
        print(f"✅ GeoJSON de estaciones exportado: {geojson_path}")
    except Exception as e_geojson_main: print(f"ERROR CRÍTICO al guardar {geojson_path}: {e_geojson_main}"); return

    try:
        gdf_unique = gdf_estaciones.drop_duplicates(subset=['code_internal'], keep='first')
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
            # Exportar solo si el código está en el archivo de estaciones
            if code in stations_processed or code not in codes_in_gdf: continue 
            
            # Seleccionar la serie, asegurarse que sea numérico y quitar NaNs
            series_to_export = pd.to_numeric(df_series[code], errors='coerce').dropna()
            
            if series_to_export.empty: continue # Saltar si no quedan datos válidos
            
            # Crear DataFrame para exportar
            df_to_export = pd.DataFrame(series_to_export)
            col_name = f"{config['display_name']}_{config['unit']}"
            df_to_export.columns = [col_name]
            
            csv_path = os.path.join(OUTPUT_DIR, f"{code}_{var_code}.csv")
            try:
                # Exportar con formato de fecha y etiqueta de índice correctos
                df_to_export.to_csv(csv_path, index=True, index_label="Fecha", date_format='%Y-%m-%d')
                count += 1
                stations_processed.add(code)
            except Exception as e_csv: print(f"ERROR al guardar CSV {csv_path}: {e_csv}")
        print(f"✅ Exportados {count} archivos CSV para '{var_code}'.")

        # Llamar a procesar_jerarquia_geoespacial CON los datos de series cargados (df_series)
        procesar_jerarquia_geoespacial(PATH_CUENCAS_SHP, ID_COLUMNA_CUENCA, 'cuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBCUENCAS_SHP, ID_COLUMNA_SUBCUENCA, 'subcuencas', gdf_estaciones, df_series, var_code, config)
        procesar_jerarquia_geoespacial(PATH_SUBSUBCUENCAS_SHP, ID_COLUMNA_SUBSUBCUENCA, 'subsubcuencas', gdf_estaciones, df_series, var_code, config)

    print("\n--- Proceso de exportación finalizado ---")

if __name__ == "__main__":
    exportar_datos_estaticos()