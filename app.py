import pandas as pd
import geopandas as gpd
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from shapely.geometry import Point
import json

# --- Configuración de FastAPI ---
app = FastAPI(
    title="API Geoespacial Hidrométrica",
    description="API que sirve datos de estaciones y cuencas hidrográficas en formato GeoJSON."
)

# Configuración de CORS
origins = ["*"] # Permite solicitudes desde cualquier origen (necesario para GitHub Pages)

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Variables Globales para Datos (Se llenan al inicio) ---
gdf_estaciones = None
gdf_cuencas = None

# --- Funciones de Carga de Datos (Se ejecutan solo al iniciar la API) ---

def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    global gdf_estaciones
    df = None
    
    # Intenta con la codificación más común
    for enc in ("cp1252", "utf-8"): 
        try:
            # **IMPORTANTE**: sep=',' y decimal='.' son correctos para tu CSV
            df = pd.read_csv("data/estaciones.csv", encoding=enc, sep=',', decimal='.')
            break
        except Exception:
            df = None
    
    if df is None:
        print("Error: No se pudo leer el archivo CSV de estaciones.")
        return None

    # Normalizar columnas a minúsculas, sin espacios, para asegurar compatibilidad
    df.columns = [c.strip().lower() for c in df.columns]

    # Asignar y convertir a numérico usando los nombres reales del CSV: 'lon' y 'lat'.
    # Usamos errors='coerce' para convertir datos inválidos (como el texto final) a NaN
    df['longitude'] = pd.to_numeric(df['lon'], errors='coerce')
    df['latitude'] = pd.to_numeric(df['lat'], errors='coerce')

    # Limpieza: Eliminar filas donde la latitud o longitud NO son números válidos (NaN)
    df_clean = df.dropna(subset=['longitude', 'latitude']).copy()

    if df_clean.empty:
        # Esto ocurre cuando todas las filas son inválidas
        print("Advertencia: No quedan estaciones con coordenadas válidas.")
        gdf_estaciones = gpd.GeoDataFrame(geometry=gpd.points_from_xy([], []), crs="EPSG:4326")
        return gdf_estaciones

    # Crear geometrías y GeoDataFrame
    geometry = [Point(xy) for xy in zip(df_clean['longitude'], df_clean['latitude'])]
    gdf_estaciones = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    print(f"[{len(gdf_estaciones)} Estaciones cargadas globalmente]")
    
    return gdf_estaciones

def cargar_cuencas():
    """Lee el Shapefile de cuencas y devuelve un GeoDataFrame."""
    global gdf_cuencas
    try:
        # Carga el Shapefile (todos los archivos auxiliares deben estar presentes)
        gdf_cuencas = gpd.read_file("cuencas_bna/cuencas_bna.shp")
        # Aseguramos que el CRS sea WGS84 para compatibilidad web
        gdf_cuencas = gdf_cuencas.to_crs("EPSG:4326")
        print(f"[{len(gdf_cuencas)} Cuencas cargadas globalmente]")
        return gdf_cuencas
    except Exception as e:
        print(f"Error al cargar cuencas.shp: {e}")
        return None

# --- EVENTO DE INICIO: Carga los datos una sola vez ---
@app.on_event("startup")
async def startup_event():
    """Ejecuta las funciones de carga al iniciar la aplicación."""
    cargar_estaciones()
    cargar_cuencas()


# --- ENDPOINT 1: ESTACIONES ---
@app.get("/api/v1/lugares", response_class=JSONResponse, tags=["Geoespacial"])
async def get_lugares():
    """Devuelve las estaciones hidrométricas en formato GeoJSON."""
    if gdf_estaciones is None or gdf_estaciones.empty:
        return JSONResponse(status_code=500, content={"error": "Datos de estaciones no disponibles"})
    
    # Conversión a GeoJSON
    estaciones_geojson = gdf_estaciones.to_json()
    # Usamos json.loads() para asegurar que se envía como objeto JSON
    return JSONResponse(content=json.loads(estaciones_geojson), media_type="application/json")


# --- ENDPOINT 2: CUENCAS ---
@app.get("/api/v1/cuencas", response_class=JSONResponse, tags=["Geoespacial"])
async def get_cuencas():
    """Devuelve las cuencas hidrográficas en formato GeoJSON."""
    if gdf_cuencas is None or gdf_cuencas.empty:
        return JSONResponse(status_code=500, content={"error": "Datos de cuencas no disponibles"})

    # Conversión a GeoJSON
    cuencas_geojson = gdf_cuencas.to_json()
    return JSONResponse(content=json.loads(cuencas_geojson), media_type="application/json")