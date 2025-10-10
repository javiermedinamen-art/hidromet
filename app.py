from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import pandas as pd
from shapely.geometry import Point
import geopandas as gpd
import json


# 1. Configuración de FastAPI
app = FastAPI(title="API Geoespacial Simple")

# Habilitar CORS para que el frontend local pueda acceder
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Permite cualquier origen (solo para desarrollo local)
    allow_methods=["*"],
    allow_headers=["*"],
)


# 2. Función de Carga y Conversión de Datos
def cargar_y_convertir_a_geojson():
    """Lee el CSV, crea geometrías y devuelve GeoJSON como un objeto Python (dict)."""
    df = None
    # Usar cp1252 o utf-8 (la solución al error anterior)
    for enc in ("cp1252", "utf-8"): 
        try:
            # Tu archivo se llama 'estaciones.csv'
            df = pd.read_csv("estaciones.csv", encoding=enc)
            break
        except Exception:
            df = None

    if df is None:
        print("Error: no se pudo leer el archivo CSV con los encodings probados.")
        return {}

    # Normalizar nombres de columnas (quitar espacios)
    df.columns = [c.strip() for c in df.columns]

    # Verificar que existan las columnas esenciales
    required_cols = ['lon', 'lat', 'name', 'code_internal', 'elevation', 'fuente']
    for col in required_cols:
        if col not in df.columns:
            print(f"Error: la columna '{col}' no se encontró en el CSV.")
            # Continuamos, pero filtramos por las que sí estén
            # En tu CSV actual, todas existen, así que esto es una medida de seguridad.
            pass 

    # Forzar a numérico, convertir strings inválidos a NaN
    df['lon'] = pd.to_numeric(df['lon'], errors='coerce')
    df['lat'] = pd.to_numeric(df['lat'], errors='coerce')

    # Eliminar filas sin coordenadas válidas
    df_clean = df.dropna(subset=['lon', 'lat']).copy()
    if df_clean.empty:
        print("Advertencia: no quedan filas con coordenadas válidas después del filtrado.")
        return {}

    # Crear geometrías y GeoDataFrame
    geometry = [Point(xy) for xy in zip(df_clean['lon'], df_clean['lat'])]
    gdf = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")

    # Eliminamos lon/lat ya que están en geometry, pero mantenemos todas las otras columnas para el popup
    try:
        gdf = gdf.drop(columns=['lon', 'lat'])
    except Exception:
        pass

    # to_json() devuelve una cadena; la cargamos como dict para enviarla como objeto JSON
    try:
        geojson_str = gdf.to_json()
        return json.loads(geojson_str)
    except Exception as e:
        print(f"Error convirtiendo GeoDataFrame a GeoJSON: {e}")
        return {}


# 3. Endpoint de la API
@app.get("/api/v1/lugares", tags=["Geoespacial"])
async def obtener_datos_geojson():
    """Devuelve los datos de lugares en formato GeoJSON (objeto JSON)."""
    geojson_data = cargar_y_convertir_a_geojson()
    if not geojson_data:
        # FastAPI y JSONResponse se encargan de la correcta serialización
        return JSONResponse(status_code=500, content={"error": "No se pudieron cargar los datos."})
    return JSONResponse(content=geojson_data)