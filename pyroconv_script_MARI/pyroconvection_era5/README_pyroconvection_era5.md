# Piroconvezione da ERA5

## 1. Installazione

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux/macOS
source .venv/bin/activate

pip install -r requirements_pyroconvection_era5.txt
```

Cartopy è opzionale. Con Conda:

```bash
conda install -c conda-forge cartopy
```

## 2. Configurazione CDS API

Registrarsi al Climate Data Store, accettare le condizioni d'uso dei due dataset ERA5 e creare il file `.cdsapirc` nella propria home:

```text
url: https://cds.climate.copernicus.eu/api
key: <PERSONAL-ACCESS-TOKEN>
```

## 3. Esecuzione di esempio

```bash
python pyroconvection_era5.py \
  --start 2026-06-20 \
  --end 2026-06-24 \
  --hours 00 03 06 09 12 15 18 21 \
  --area 45.2 8.8 41.4 13.5 \
  --outdir output_piroconvezione
```

Su PowerShell usare il backtick `` ` `` al posto di `\` per andare a capo, oppure scrivere il comando su una sola riga.

## 4. Output

- `pyroconvection_era5.nc`: dataset completo con classe e diagnostiche.
- `geotiff/pyro_class_*.tif`: raster EPSG:4326 per ogni ora.
- `png/pyro_class_*.png`: mappe singole.
- `pyroconvection_panel.png`: tavola multi-pannello.
- `classi_piroconvezione.csv`: codici, etichette e colori.

Classi:

- 0: nessuna classe / condizioni insufficienti
- 1: Convective plumes
- 2: Overshooting PyroCu potenziale
- 3: Resilient PyroCu potenziale
- 4: PyroCb / deep PyroCu potenziale

## 5. Nota scientifica

ERA5 è una rianalisi. Lo script produce una diagnosi della predisposizione atmosferica, non una previsione operativa né una probabilità di occorrenza. Le classi 2-4 presuppongono un incendio sufficientemente intenso. Le soglie sono iniziali e devono essere validate su casi locali.

## 6. Confini regionali nei PNG

Per aggiungere i confini amministrativi regionali alle mappe singole e alla
tavola multi-pannello usare:

```powershell
python pyroconvection_era5.py `
  --start 2026-06-20 `
  --end 2026-06-24 `
  --hours 00 03 06 09 12 15 18 21 `
  --area 45.2 8.8 41.4 13.5 `
  --regional-boundaries `
  --outdir output_piroconvezione
```

Lo script usa Natural Earth `admin_1_states_provinces_lines` a scala 1:10m.
Al primo utilizzo Cartopy può scaricare il layer nella propria cache.

È possibile utilizzare un proprio shapefile geografico EPSG:4326:

```powershell
python pyroconvection_era5.py ... `
  --regions-shapefile "C:\dati\confini_regioni_wgs84.shp"
```

I confini sono sovrapposti soltanto ai PNG; i GeoTIFF restano raster tematici
senza linee vettoriali incorporate.
