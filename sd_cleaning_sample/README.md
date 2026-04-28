# sd_cleaning (muestra: baseline)

Código **mínimo** que aísla **un solo paso** del pipeline de profundidad de nieve: la corrección de baseline en segmentos de baja varianza (`detect_flat_offsets` en `sd_ground_baseline.py`). El resto del procesamiento está fuera de esta carpeta.

## Contenido

| Archivo | Rol |
|---------|-----|
| `sd_ground_baseline.py` | Algoritmo + CLI opcional para visualizar |
| `baseline_viz.ipynb` | Misma idea paso a paso en Jupyter |
| `SD_ValleOlivares_05706003.csv` | Ejemplo reproducible (`date` + serie `05706003`) |
| `requirements.txt` | Dependencias |

## Uso rápido

```bash
python -m pip install -r requirements.txt
```

- **Notebook**: abrir `baseline_viz.ipynb` y ejecutar las celdas en orden (ajusta la ruta del CSV si hace falta).
- **Línea de comandos**: con el CSV de ejemplo por defecto:

```bash
python sd_ground_baseline.py
```

Otros CSV necesitan una columna `date` y la columna de valores (única numérica salvo fecha, o bien `--col`). Sin entorno gráfico, puedes guardar la figura con `-o salida.png`.
