import pandas as pd
import geopandas as gpd
from fastapi import FastAPI
from fastapi.responses import JSONResponse, HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from shapely.geometry import Point
import json
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from datetime import datetime

# --- Configuración de FastAPI ---
app = FastAPI(
    title="API Geoespacial Hidrométrica (Estaciones)",
    description="API que sirve datos de estaciones hidrométricas en formato GeoJSON y genera series temporales con Plotly."
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
gdf_estaciones = None 
df_series = None 

# --- Funciones de Carga de Datos (Sin Cambios Relevantes) ---

def cargar_estaciones():
    """Lee el CSV de estaciones y devuelve un GeoDataFrame."""
    global gdf_estaciones
    # ... (código de cargar_estaciones permanece igual)
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
    
    required_cols = ['lon', 'lat', 'code_internal', 'name', 'fuente', 'basin']
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
    gdf_estaciones = gpd.GeoDataFrame(df_clean, geometry=geometry, crs="EPSG:4326")
    print(f"[{len(gdf_estaciones)} Estaciones cargadas globalmente]")
    
    return gdf_estaciones

# --- Función de Carga de Datos (Actualizada con formato forzado) ---

def cargar_series_temporales():
    """Carga los datos de series temporales, asegura la conversión a datetime y numérico."""
    global df_series
    file_path = "data/PPT_clean_v1.csv"
    
    print("--- INICIANDO CARGA DE SERIES TEMPORALES ---")
    
    try:
        # 1. Lectura inicial con separadores y decimales forzados (asumiendo ',' y '.')
        df = pd.read_csv(
            file_path, 
            encoding="utf-8", 
            sep=',',              
            decimal='.'          
        )
        
        # 2. Limpieza y conversión a minúsculas para asegurar coincidencia
        df.columns = df.columns.astype(str).str.strip() 
        df.rename(columns={c: c.lower() for c in df.columns}, inplace=True)
        
        # DEPURA: ¿Qué columnas ve Pandas justo después de leer?
        print(f"DEPURACIÓN: Columnas leídas y limpiadas: {df.columns.tolist()}")

        # 3. Convertir la columna 'date' al índice de tiempo
        if 'date' not in df.columns:
            print("Error CRÍTICO: La columna 'date' no se encontró en el CSV.")
            return False

        # Forzar conversión de fecha, eliminando filas con fechas no válidas
        df['date'] = pd.to_datetime(df['date'].astype(str).str.strip(), errors='coerce')
        df.dropna(subset=['date'], inplace=True) 

        # 4. Establecer el índice y convertir todas las columnas de datos a numérico
        df_series = df.set_index('date')
        
        for col in df_series.columns:
            # Forzamos la conversión a numérico (deja NaN si no puede convertir)
            df_series[col] = pd.to_numeric(df_series[col], errors='coerce')
        
        # 5. Limpieza final: Series que quedaron solo con NaN se eliminan
        df_series.dropna(axis=1, how='all', inplace=True) 
        
        if df_series.empty or len(df_series.columns) == 0:
            print("Advertencia CRÍTICA: La serie temporal está vacía. Revise los datos y encabezados.")
            return False

        print(f"[{len(df_series.columns)} Series Temporales cargadas (¡Éxito!)]")
        # DEPURA: Muestra las primeras filas y los tipos de datos
        print("\nDEPURACIÓN: Primeras filas del DataFrame final (df_series.head()):\n", df_series.head())
        print("\nDEPURACIÓN: Tipos de datos de columnas:\n", df_series.dtypes)
        return True
        
    except FileNotFoundError:
        print(f"Error FATAL: Archivo no encontrado en '{file_path}'.")
        return False
    except Exception as e:
        print(f"Error FATAL INESPERADO al procesar '{file_path}': {e}")
        return False

# --- EVENTO DE INICIO ---
@app.on_event("startup")
async def startup_event():
    cargar_estaciones()
    cargar_series_temporales() 

# --- ENDPOINTS ---

@app.get("/api/v1/lugares", response_class=JSONResponse, tags=["Geoespacial"])
async def get_lugares():
    """Devuelve las estaciones hidrométricas en formato GeoJSON."""
    if gdf_estaciones is None or gdf_estaciones.empty:
        return JSONResponse(status_code=500, content={"error": "Datos de estaciones no disponibles"})
    
    estaciones_geojson = gdf_estaciones.to_json()
    return JSONResponse(content=json.loads(estaciones_geojson), media_type="application/json")


# 🟢 ENDPOINT CON GENERACIÓN DE BOTONES DE ZOOM
@app.get("/api/v1/series/chart/{code_internal}", response_class=HTMLResponse, tags=["Datos"])
async def get_series_chart(code_internal: str):
    """Genera y devuelve el gráfico interactivo (HTML) con botones de zoom anual."""
    global df_series, gdf_estaciones
    
    if df_series is None or df_series.empty:
        return HTMLResponse(content="<div>Error: Datos de series temporales no disponibles.</div>", status_code=500)
    
    code_internal = str(code_internal).strip()

    if code_internal not in df_series.columns:
        return HTMLResponse(content=f"<div>Error: Serie temporal no encontrada para el código: {code_internal}</div>", status_code=404)
    
    # 1. Preparar datos y nombre
    series_data = df_series[[code_internal]].dropna()
    
    station_info = gdf_estaciones[gdf_estaciones['code_internal'] == code_internal]
    station_name = station_info['name'].iloc[0] if not station_info.empty else code_internal
    
    # 2. Crear el gráfico base (Plotly Graph Objects para control avanzado)
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=series_data.index,
        y=series_data[code_internal],
        mode='lines',
        name='Precipitación',
        line=dict(color='#007bff', width=2)
    ))

    # 3. Generar la lista de años únicos
    # Usamos el índice de la serie temporal (que es la fecha)
    unique_years = sorted(series_data.index.year.unique().tolist())
    
    updatemenus = [
        # Menú principal de zoom
        dict(
            type="buttons",
            direction="right",
            active=0, # Por defecto, la serie completa es el estado activo (zoom 0)
            x=0.0, 
            y=1.1, 
            xanchor="left",
            yanchor="top",
            buttons=[
                # Botón de Serie Completa (Reset Zoom)
                dict(
                    label="Serie Completa",
                    method="relayout",
                    args=[{"xaxis.range": [None, None]}] # Establece el rango a automático
                )
            ] + [
                # Botones por cada año
                dict(
                    label=str(year),
                    method="relayout",
                    args=[{
                        # Establece el rango de fecha desde el 1 de Enero al 31 de Diciembre del año
                        "xaxis.range": [
                            datetime(year, 1, 1).strftime('%Y-%m-%d'), 
                            datetime(year, 12, 31).strftime('%Y-%m-%d')
                        ]
                    }]
                ) for year in unique_years
            ]
        )
    ]
    
    # 4. Configurar Layout
    fig.update_layout(
        title={
            'text': f'Precipitación diaria - {station_name}',
            'y':0.95,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        updatemenus=updatemenus,
        margin=dict(l=40, r=40, t=100, b=40),
        xaxis_title="Fecha",
        yaxis_title="Precipitación (mm)",
        hovermode="x unified",
        xaxis=dict(
            rangeselector=dict(
                buttons=list([
                    dict(count=1, label="1m", step="month", stepmode="backward"),
                    dict(count=6, label="6m", step="month", stepmode="backward"),
                    dict(count=1, label="YTD", step="year", stepmode="todate"),
                    dict(count=1, label="1a", step="year", stepmode="backward"),
                    dict(step="all")
                ])
            ),
            rangeslider=dict(visible=True),
            type="date"
        )
    )

    # Convertir el gráfico a HTML interactivo
    chart_html = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')
    
    return HTMLResponse(content=chart_html)