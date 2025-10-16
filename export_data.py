import pandas as pd
import geopandas as gpd
import os
import json
import sys
from shapely.geometry import Point

# --------------------------------------------------------------------------
# --- CONFIGURACIÓN DE RUTAS Y COLUMNAS ---
# ¡¡AJUSTA ESTAS RUTAS Y NOMBRES DE COLUMNAS SEGÚN TUS ARCHIVOS!!
# --------------------------------------------------------------------------
OUTPUT_DIR = "data_static"

# Ruta al shapefile de Cuencas y nombre de su columna de ID
PATH_CUENCAS_SHP = "data/cuencas_bna/cuencas_bna.shp"
ID_COLUMNA_CUENCA = 'COD_CUEN'
NOMBRE_COLUMNA_CUENCA = 'NOM_CUEN'

# Ruta al shapefile de Subcuencas y nombre de su columna de ID
PATH_SUBCUENCAS_SHP = "data/subcuencas_bna/subcuencas_bna.shp"
ID_COLUMNA_SUBCUENCA = 'COD_SUBC'
NOMBRE_COLUMNA_SUBCUENCA = 'NOM_SUBC'

# Ruta al shapefile de Subsubcuencas y nombres de sus columnas
PATH_SUBSUBCUENCAS_SHP = "data/subsubcuencas_bna/subsubcuencas_bna.shp" 
ID_COLUMNA_SUBSUBCUENCA = 'COD_SSUBC' 
NOMBRE_COLUMNA_SUBSUBCUENCA = 'NOM_SSUBC'

# --------------------------------------------------------------------------

def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    possible_paths = [
        "data/estaciones.csv",
        "data\\estaciones.csv",
        "estaciones.csv",
    ]

    df = None
    last_exc = None
    for p in possible_paths:
        if os.path.exists(p):
            try:
                df = pd.read_csv(p, encoding="utf-8", sep=',', decimal='.')
                print(f"Leyendo estaciones desde: {p}")
                break
            except Exception:
                try:
                    df = pd.read_csv(p, encoding="cp1252", sep=',', decimal='.')
                    print(f"Leyendo estaciones desde: {p} con cp1252")
                    break
                except Exception as e2:
                    last_exc = e2
    if df is None:
        print(f"Error CRÍTICO al leer estaciones.csv. Intentados: {possible_paths}. Último error: {last_exc}")
        return None

    df.columns = [c.strip().lower() for c in df.columns]
    required_cols = ['lon', 'lat', 'code_internal', 'name', 'fuente', 'basin']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Faltan columnas requeridas en estaciones.csv. Necesarias: {required_cols}")
        return None

    df['longitude'] = pd.to_numeric(df['lon'].astype(str).str.strip(), errors='coerce')
    df['latitude'] = pd.to_numeric(df['lat'].astype(str).str.strip(), errors='coerce')
    df_clean = df.dropna(subset=['longitude', 'latitude']).copy()

    if df_clean.empty:
        print("Advertencia: No quedan estaciones con coordenadas válidas.")
        return gpd.GeoDataFrame()

    geometry = [Point(xy) for xy in zip(df_clean['longitude'], df_clean['latitude'])]
    gdf = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    return gdf

def cargar_series_temporales():
    """Carga los datos de series temporales, asegura la conversión a datetime y numérico."""
    possible_paths = [
        "data/PPT_clean_v1.csv",
        "data\\PPT_clean_v1.csv",
        "PPT_clean_v1.csv",
    ]
    file_path = None
    for p in possible_paths:
        if os.path.exists(p):
            file_path = p
            break
    if file_path is None:
        print(f"Error: no se encontró el archivo de series. Buscados: {possible_paths}")
        return None

    try:
        df = pd.read_csv(file_path, encoding="utf-8", sep=',', decimal='.')
        df.columns = df.columns.astype(str).str.strip()
        df.rename(columns={c: c.lower() for c in df.columns}, inplace=True)

        if 'date' not in df.columns:
            print("Error CRÍTICO: La columna 'date' no se encontró en el CSV.")
            return None

        df['date'] = pd.to_datetime(df['date'].astype(str).str.strip(), errors='coerce')
        df.dropna(subset=['date'], inplace=True)

        df_series = df.set_index('date')

        for col in df_series.columns:
            df_series[col] = pd.to_numeric(df_series[col], errors='coerce')

        df_series.dropna(axis=1, how='all', inplace=True)
        return df_series

    except Exception as e:
        print(f"Error al cargar series temporales: {e}")
        return None

def procesar_jerarquia_geoespacial(shapefile_path, id_columna, output_prefix, gdf_estaciones, df_series):
    """
    Procesa una capa geoespacial.
    Añade la cuenta de estaciones a cada polígono.
    """
    print(f"\n--- Procesando jerarquía: '{output_prefix}' ---")

    if not os.path.exists(shapefile_path):
        print(f"Advertencia: No se encontró el shapefile en '{shapefile_path}'. Saltando esta jerarquía.")
        return

    try:
        gdf_jerarquia = gpd.read_file(shapefile_path)
        gdf_jerarquia = gdf_jerarquia.to_crs(gdf_estaciones.crs)
        if id_columna not in gdf_jerarquia.columns:
            print(f"Error CRÍTICO: La columna '{id_columna}' no existe en '{shapefile_path}'.")
            return
    except Exception as e:
        print(f"Error CRÍTICO al leer '{shapefile_path}': {e}")
        return

    gdf_estaciones_con_jerarquia = gpd.sjoin(gdf_estaciones, gdf_jerarquia, how="inner", predicate="within")
    ids_con_estaciones = gdf_estaciones_con_jerarquia[id_columna].unique()
    gdf_filtrado = gdf_jerarquia[gdf_jerarquia[id_columna].isin(ids_con_estaciones)].copy()

    if gdf_filtrado.empty:
        print(f"Advertencia: Ninguna estación encontrada para la jerarquía '{output_prefix}'.")
        return

    # Contar estaciones por polígono y añadirlo al GeoDataFrame
    station_counts = gdf_estaciones_con_jerarquia.groupby(id_columna).size()
    station_counts.name = 'n_estaciones'
    gdf_filtrado = gdf_filtrado.merge(station_counts, left_on=id_columna, right_index=True, how='left')
    gdf_filtrado['n_estaciones'] = gdf_filtrado['n_estaciones'].fillna(0).astype(int)
    print(f"✅ Se calculó el número de estaciones para {len(gdf_filtrado)} polígonos.")

    geojson_path = os.path.join(OUTPUT_DIR, f"{output_prefix}.geojson")
    gdf_filtrado.to_file(geojson_path, driver="GeoJSON")
    print(f"✅ GeoJSON de '{output_prefix}' con cuenta de estaciones exportado a: {geojson_path}")

    count_poligonos = 0
    for poly_id in ids_con_estaciones:
        estaciones_en_poligono = gdf_estaciones_con_jerarquia[gdf_estaciones_con_jerarquia[id_columna] == poly_id].copy()
        
        estaciones_en_poligono['elevation'] = pd.to_numeric(estaciones_en_poligono['elevation'], errors='coerce')
        
        estaciones_ordenadas = estaciones_en_poligono.sort_values(by='elevation', ascending=False, na_position='last')
        
        codigos_estaciones_ordenados = estaciones_ordenadas['code_internal'].tolist()
        df_series_poligono = df_series.reindex(columns=codigos_estaciones_ordenados)
        
        monthly_agg = df_series_poligono.resample('MS').agg(['sum', 'count'])

        for col in codigos_estaciones_ordenados:
            if (col, 'count') in monthly_agg.columns:
                days_in_month = monthly_agg.index.days_in_month
                condition = monthly_agg[(col, 'count')] < days_in_month * 0.8
                monthly_agg.loc[condition, (col, 'sum')] = None

        df_sum_monthly = monthly_agg.xs('sum', level=1, axis=1).dropna(how='all')

        if df_sum_monthly.empty:
            continue
        
        start_date = df_sum_monthly.index.min()
        end_date = df_sum_monthly.index.max()
        
        complete_date_index = pd.date_range(start_date, end_date, freq='MS')
        df_complete = df_sum_monthly.reindex(complete_date_index)

        output_json = {}
        for station_code in codigos_estaciones_ordenados:
            if station_code not in df_complete.columns:
                continue

            station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
            
            values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_complete[station_code]]
            elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']
            
            output_json[station_code] = {
                'name': station_info['name'],
                'elevation': elevation_value,
                'dates': df_complete.index.strftime('%Y-%m-%d').tolist(),
                'values': values_with_nulls
            }

        json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}.json")
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(output_json, f, ensure_ascii=False)

        # --- INICIO DEL NUEVO CÓDIGO PARA CÁLCULO ANUAL ROBUSTO ---
        # 1. Calcular el perfil climatológico y los Meses Hidrológicamente Significativos (MHS) para cada estación
        hsm_por_estacion = {}
        for station_code in df_sum_monthly.columns:
            # Agrupar por mes y calcular la media a largo plazo para esta estación
            climatologia_mensual = df_sum_monthly[station_code].groupby(df_sum_monthly.index.month).mean()
            
            # Ordenar los meses por su contribución de mayor a menor
            climatologia_ordenada = climatologia_mensual.sort_values(ascending=False)
            
            # Calcular el umbral del 75% de la precipitación total anual promedio
            umbral_precipitacion = climatologia_mensual.sum() * 0.75
            
            # Encontrar los meses que superan este umbral acumulado
            suma_acumulada = climatologia_ordenada.cumsum()
            meses_significativos = suma_acumulada[suma_acumulada <= umbral_precipitacion].index.tolist()
            
            # Aseguramos que al menos el mes más lluvioso esté, incluso si por sí solo supera el 75%
            if not meses_significativos and not climatologia_ordenada.empty:
                meses_significativos = [climatologia_ordenada.index[0]]
                
            hsm_por_estacion[station_code] = meses_significativos

        # 2. Verificar la validez de cada año para cada estación según las nuevas reglas
        df_annual = df_sum_monthly.resample('YS').sum(min_count=1) # Suma anual preliminar
        monthly_counts = df_sum_monthly.resample('YS').count() # Conteo de meses válidos por año

        for station_code in df_annual.columns:
            # Iterar sobre cada año en el índice del DataFrame anual
            for year_start_date in df_annual.index:
                year = year_start_date.year
                
                # Condición A: Verificar si el año tiene al menos 9 meses válidos
                try:
                    num_meses_validos = monthly_counts.loc[year_start_date, station_code]
                except KeyError:
                    num_meses_validos = 0
                    
                if num_meses_validos < 9:
                    df_annual.loc[year_start_date, station_code] = None # Anular si no cumple
                    continue # Pasar al siguiente año para esta estación

                # Condición B: Verificar si TODOS los MHS están presentes en ese año
                meses_clave = hsm_por_estacion.get(station_code, [])
                if not meses_clave: # Si no hay meses clave, continuar
                    continue
                
                # Extraer los datos mensuales solo para el año actual y la estación actual
                datos_del_ano = df_sum_monthly.loc[str(year), station_code]
                
                # Verificar si alguno de los meses clave es NULO en ese año
                for mes_clave in meses_clave:
                    if mes_clave not in datos_del_ano.index.month or pd.isna(datos_del_ano[datos_del_ano.index.month == mes_clave].iloc[0]):
                        df_annual.loc[year_start_date, station_code] = None # Anular si falta un mes clave
                        break # No es necesario seguir revisando otros meses clave para este año

        if not df_annual.empty:
            output_annual_json = {}
            for station_code in codigos_estaciones_ordenados:
                if station_code not in df_annual.columns:
                    continue

                station_info = estaciones_ordenadas.set_index('code_internal').loc[station_code]
                
                values_with_nulls = [None if pd.isna(val) else round(val, 1) for val in df_annual[station_code]]
                elevation_value = None if pd.isna(station_info['elevation']) else station_info['elevation']

                output_annual_json[station_code] = {
                    'name': station_info['name'],
                    'elevation': elevation_value,
                    'years': df_annual.index.year.tolist(),
                    'values': values_with_nulls
                }
            
            annual_json_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_{poly_id}_annual.json")
            with open(annual_json_path, 'w', encoding='utf-8') as f:
                json.dump(output_annual_json, f, ensure_ascii=False)
        # --- FIN DEL NUEVO CÓDIGO ---
        
        count_poligonos += 1

    print(f"✅ Exportados {count_poligonos} archivos JSON de datos para '{output_prefix}'.")
    
def exportar_datos_estaticos():
    """
    Función principal que orquesta la carga de datos y la exportación de todos
    los archivos estáticos necesarios para el visualizador.
    """
    print("Iniciando exportación de datos estáticos...")
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # --- PASO 1: Cargar datos base ---
    gdf_estaciones = cargar_estaciones()
    if gdf_estaciones is None or gdf_estaciones.empty:
        print("Fallo la carga de estaciones. Abortando.")
        return

    df_series = cargar_series_temporales()
    if df_series is None:
        print("Fallo la carga de series temporales. Abortando.")
        return

    # --- PASO 2: Exportar datos de estaciones individuales ---
    print("\n--- Exportando datos de estaciones individuales ---")
    geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    gdf_estaciones.to_file(geojson_path, driver="GeoJSON")
    print(f"✅ GeoJSON de estaciones exportado a: {geojson_path}")

    count = 0
    for code in df_series.columns:
        series_to_export = df_series[[code]].copy().dropna()
        series_to_export.columns = ['Precipitacion_mm']
        csv_path = os.path.join(OUTPUT_DIR, f"{code}.csv")
        series_to_export.to_csv(csv_path, index=True, index_label="Fecha")
        count += 1
    print(f"✅ Exportados {count} archivos CSV de series temporales.")

    try:
        station_names = gdf_estaciones.set_index('code_internal')['name'].to_dict()
        with open(os.path.join(OUTPUT_DIR, "station_names.json"), 'w', encoding='utf-8') as f:
            json.dump(station_names, f, ensure_ascii=False, indent=2)
        print("✅ Archivo de nombres de estaciones exportado.")
    except Exception as e:
        print(f"Advertencia: no se pudo exportar station_names.json: {e}")

    # --- PASO 3: Procesar y exportar datos por jerarquía geoespacial ---
    # 3.1 Procesar Cuencas
    procesar_jerarquia_geoespacial(
        shapefile_path=PATH_CUENCAS_SHP,
        id_columna=ID_COLUMNA_CUENCA,
        output_prefix='cuencas',
        gdf_estaciones=gdf_estaciones,
        df_series=df_series
    )

    # 3.2 Procesar Subcuencas
    procesar_jerarquia_geoespacial(
        shapefile_path=PATH_SUBCUENCAS_SHP,
        id_columna=ID_COLUMNA_SUBCUENCA,
        output_prefix='subcuencas',
        gdf_estaciones=gdf_estaciones,
        df_series=df_series
    )

    # 3.3 Procesar Subsubcuencas
    procesar_jerarquia_geoespacial(
        shapefile_path=PATH_SUBSUBCUENCAS_SHP,
        id_columna=ID_COLUMNA_SUBSUBCUENCA,
        output_prefix='subsubcuencas',
        gdf_estaciones=gdf_estaciones,
        df_series=df_series
    )

    print("\n--- Proceso de exportación finalizado ---")


if __name__ == "__main__":
    exportar_datos_estaticos()