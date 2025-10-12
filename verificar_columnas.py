import geopandas as gpd

# --- ¡¡AJUSTA ESTA RUTA!! ---
ruta_a_tu_shapefile = "data/cuencas_bna/cuencas_bna.shp"

try:
    gdf = gpd.read_file(ruta_a_tu_shapefile)
    print(f"Las columnas disponibles en '{ruta_a_tu_shapefile}' son:")
    # Imprimimos la lista de columnas para que puedas ver el nombre exacto
    print(list(gdf.columns))
except Exception as e:
    print(f"No se pudo leer el archivo: {e}")