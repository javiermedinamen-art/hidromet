import pandas as pd
import geopandas as gpd
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from shapely.geometry import Point
import json

# --- Configuración de FastAPI ---
app = FastAPI(
    title="API Geoespacial Hidrométrica (Estaciones)",
    description="API que sirve datos de estaciones hidrométricas en formato GeoJSON."
)

origins = ["*"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Variables Globales para Datos ---
# Solo necesitamos gdf_estaciones
gdf_estaciones = None 

# --- Funciones de Carga de Datos ---

def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    global gdf_estaciones
    df = None
    
    try:
        # Intenta leer con UTF-8
        df = pd.read_csv("data/estaciones.csv", encoding="utf-8", sep=',', decimal='.')
    except Exception:
        try:
            # Intenta leer con CP1252 si falla UTF-8
            df = pd.read_csv("data/estaciones.csv", encoding="cp1252", sep=',', decimal='.')
        except Exception as e:
            print(f"Error CRÍTICO al leer estaciones.csv: {e}")
            return None

    # Normalización de nombres de columna: minúsculas y sin espacios
    df.columns = [c.strip().lower() for c in df.columns]
    
    # Asegurarse de que las columnas críticas existan
    required_cols = ['lon', 'lat', 'nom_cuen', 'fuente']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Faltan columnas requeridas en estaciones.csv. Necesarias: {required_cols}")
        return None

    # Conversión numérica de coordenadas
    df['longitude'] = pd.to_numeric(df['lon'].astype(str).str.strip(), errors='coerce')
    df['latitude'] = pd.to_numeric(df['lat'].astype(str).str.strip(), errors='coerce')

    df_clean = df.dropna(subset=['longitude', 'latitude']).copy()
    
    if df_clean.empty:
        print("Advertencia: No quedan estaciones con coordenadas válidas.")
        return gpd.GeoDataFrame(geometry=gpd.points_from_xy([], []), crs="EPSG:4326")

    # Crear GeoDataFrame
    geometry = [Point(xy) for xy in zip(df_clean['longitude'], df_clean['latitude'])]
    gdf_estaciones = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    print(f"[{len(gdf_estaciones)} Estaciones cargadas globalmente]")
    
    return gdf_estaciones

# --- EVENTO DE INICIO ---
@app.on_event("startup")
async def startup_event():
    cargar_estaciones()

# --- ENDPOINTS ---

@app.get("/api/v1/lugares", response_class=JSONResponse, tags=["Geoespacial"])
async def get_lugares():
    """Devuelve las estaciones hidrométricas en formato GeoJSON."""
    if gdf_estaciones is None or gdf_estaciones.empty:
        return JSONResponse(status_code=500, content={"error": "Datos de estaciones no disponibles"})
    
    # Conversión a GeoJSON
    estaciones_geojson = gdf_estaciones.to_json()
    return JSONResponse(content=json.loads(estaciones_geojson), media_type="application/json")

# Se eliminan los endpoints /api/v1/cuencas y /api/v1/cuencas/filtro