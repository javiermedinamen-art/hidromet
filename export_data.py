import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import os
import json

# --- Configuración de Rutas de Salida ---
# Los archivos estáticos se guardarán en una carpeta 'data_static'
OUTPUT_DIR = "data_static"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# --- Funciones de Carga de Datos (Copiadas desde app.py) ---
# (Asume que estas funciones usan tus archivos 'data/estaciones.csv' y 'data/PPT_clean_v1.csv')

def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    # ... (Cuerpo de la función cargar_estaciones copiado de tu app.py)
    df = None
    
    try:
        df = pd.read_csv("data/estaciones.csv", encoding="utf-8", sep=',', decimal='.')
    except Exception:
        try:
            df = pd.read_csv("data/estaciones.csv", encoding="cp1252", sep=',', decimal='.')
        except Exception as e:
            print(f"Error CRÍTICO al leer estaciones.csv: {e}")
            return None

    df.columns = [c.strip().lower() for c in df.columns]
    required_cols = ['lon', 'lat', 'code_internal', 'name', 'fuente', 'basin', 'elevation']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Faltan columnas requeridas en estaciones.csv. Necesarias: {required_cols}")
        return None

    df['longitude'] = pd.to_numeric(df['lon'].astype(str).str.strip(), errors='coerce')
    df['latitude'] = pd.to_numeric(df['lat'].astype(str).str.strip(), errors='coerce')
    df_clean = df.dropna(subset=['longitude', 'latitude']).copy()
    
    if df_clean.empty:
        print("Advertencia: No quedan estaciones con coordenadas válidas.")
        return gpd.GeoDataFrame(geometry=gpd.points_from_xy([], []), crs="EPSG:4326")

    geometry = [Point(xy) for xy in zip(df_clean['longitude'], df_clean['latitude'])]
    gdf = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    return gdf

def cargar_series_temporales(): 
    """Carga los datos de series temporales, asegura la conversión a datetime y numérico."""
    # ... (Cuerpo de la función cargar_series_temporales copiado de tu app.py)
    file_path = "data/PPT_clean_v1.csv"
    
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

# --- Función Principal de Exportación ---

def exportar_datos_estaticos():
    print("Iniciando exportación de datos estáticos...")
    
    # 1. Cargar Estaciones y exportar GeoJSON principal
    gdf_estaciones = cargar_estaciones()
    if gdf_estaciones is None or gdf_estaciones.empty:
        print("Fallo la carga de estaciones. Abortando.")
        return

    # Exportar el GeoJSON principal que usa Leaflet para dibujar
    geojson_path = os.path.join(OUTPUT_DIR, "estaciones.geojson")
    gdf_estaciones.to_file(geojson_path, driver="GeoJSON")
    print(f"✅ 1. GeoJSON principal exportado a: {geojson_path}")
    
    # 2. Cargar Series Temporales
    df_series = cargar_series_temporales()
    if df_series is None:
        print("Fallo la carga de series temporales. Abortando exportación de series.")
        return

    # 3. Exportar un CSV individual por cada estación
    count = 0
    for code in df_series.columns:
        # Renombrar la columna a 'Precipitacion' antes de exportar
        series_to_export = df_series[[code]].copy().dropna()
        series_to_export.columns = ['Precipitacion_mm'] 
        
        csv_path = os.path.join(OUTPUT_DIR, f"{code}.csv")
        # Exportamos con el índice (la fecha)
        series_to_export.to_csv(csv_path, index=True, index_label="Fecha") 
        count += 1

    print(f"✅ 2. Exportados {count} archivos CSV de series temporales a la carpeta: {OUTPUT_DIR}")
    
    # Opcional: Crear un archivo de mapeo para los nombres
    station_names = gdf_estaciones.set_index('code_internal')['name'].to_dict()
    with open(os.path.join(OUTPUT_DIR, "station_names.json"), 'w') as f:
        json.dump(station_names, f)
    print(f"✅ 3. Archivo de nombres de estaciones exportado a: {os.path.join(OUTPUT_DIR, 'station_names.json')}")

if __name__ == "__main__":
    exportar_datos_estaticos()