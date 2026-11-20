# Piroconvezione da ICON-2I / MeteoHub

Questo pacchetto adatta la diagnostica originariamente sviluppata per ERA5 ai
campi previsionali **ICON-2I** pubblicati da MeteoHub/ItaliaMeteo.

## Modalità implementata

Lo script usa il repository open-data pubblico:

`https://meteohub.agenziaitaliameteo.it/nwp/ICON-2I_SURFACE_PRESSURE_LEVELS/`

Non sono necessarie credenziali. Il programma:

1. legge l'indice dei run disponibili;
2. sceglie il run più recente completo, oppure quello richiesto;
3. scarica i GRIB delle variabili necessarie;
4. apre i GRIB con `cfgrib`/ecCodes;
5. seleziona scadenze, livelli e area;
6. calcola la diagnostica di piroconvezione;
7. produce NetCDF, GeoTIFF e PNG.

## Variabili ICON-2I utilizzate

Superficie:

- `T_2M`
- `TD_2M`
- `U_10M`
- `V_10M`
- `PS`
- `HSURF`
- `FR_LAND`

Livelli isobarici:

- `FI`
- `T`
- `QV`
- `U`
- `V`

`FI` viene trattato come geopotenziale e convertito in quota dividendo per `g`.
`HSURF` viene trattato come altezza geometrica del terreno.

ICON-2I non fornisce nel repository open-data usato dallo script un campo BLH
equivalente a quello ERA5. L'altezza dell'ABL viene quindi stimata dal profilo
mediante Bulk Richardson number.

## Installazione consigliata

Su Windows è preferibile usare conda-forge, perché installa anche ecCodes:

```powershell
conda create -n pyro_icon2i -c conda-forge `
  python=3.11 numpy pandas xarray netcdf4 requests `
  cfgrib eccodes rasterio matplotlib cartopy

conda activate pyro_icon2i
```

In alternativa:

```bash
pip install -r requirements_pyroconvection_icon2i.txt
```

Con `pip`, verificare che ecCodes sia installato e correttamente visibile a
`cfgrib`.

## Elenco dei run disponibili

```powershell
python pyroconvection_icon2i_meteohub.py --list-runs
```

## Esecuzione sul dominio Toscana

Primo test consigliato, scadenze ogni 3 ore fino a +48 h e griglia
sottocampionata di un fattore 2:

```powershell
python pyroconvection_icon2i_meteohub.py `
  --run latest `
  --lead-hours 0 3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 `
  --area 45.2 8.8 41.4 13.5 `
  --stride 2 `
  --outdir output_piroconvezione_icon2i
```

Alla risoluzione nativa:

```powershell
python pyroconvection_icon2i_meteohub.py `
  --run latest `
  --lead-hours 0 3 6 9 12 15 18 21 24 `
  --area 45.2 8.8 41.4 13.5 `
  --stride 1 `
  --outdir output_piroconvezione_icon2i
```

## Riutilizzo dei file già scaricati

```powershell
python pyroconvection_icon2i_meteohub.py `
  --no-download `
  --run 2026071112 `
  --lead-hours 0 3 6 9 12 15 18 21 24 `
  --area 45.2 8.8 41.4 13.5 `
  --outdir output_piroconvezione_icon2i
```

I GRIB sono attesi in:

```text
OUTDIR/
└── raw/
    └── YYYYMMDDHH/
        ├── T_2M/
        ├── TD_2M/
        ├── U_10M/
        ├── V_10M/
        ├── PS/
        ├── HSURF/
        ├── FR_LAND/
        ├── FI/
        ├── T/
        ├── QV/
        ├── U/
        └── V/
```

## Output

```text
OUTDIR/
├── pyroconvection_icon2i_YYYYMMDDHH.nc
├── pyroconvection_icon2i_panel_YYYYMMDDHH.png
├── classi_piroconvezione.csv
├── geotiff/
│   └── pyro_class_YYYYMMDDTHHMMZ.tif
└── png/
    └── pyro_class_YYYYMMDDTHHMMZ.png
```

## Nota su dimensione e tempi

La modalità open-data scarica i file GRIB completi del dominio ICON-2I. I
campi ai livelli isobarici possono essere molto pesanti. L'opzione `--stride`
riduce il costo dell'elaborazione, ma non il volume scaricato.

Per una pipeline operativa più efficiente conviene usare l'API autenticata di
estrazione MeteoHub, applicando direttamente filtri su area, livelli,
variabili e scadenze. Per implementare tale modalità serve il JSON di una
richiesta ICON-2I creata nell'interfaccia MeteoHub tramite:

`Dataset > Estrazione dati > Le mie richieste > Copy to clipboard`

Non inserire username o password nel JSON né inviarli insieme al file.

## Limiti scientifici

Le classi 2-4 presuppongono un incendio sufficientemente intenso. Il prodotto
rappresenta predisposizione atmosferica, non probabilità calibrata di
formazione di PyroCu/PyroCb. Le soglie devono essere validate sui casi locali.


## Correzione v1.1.0 — livelli isobarici MeteoHub

Nei GRIB ICON-2I di MeteoHub ogni file ai livelli di pressione contiene un
singolo livello. `cfgrib` lo restituisce normalmente come coordinata scalare
`isobaricInhPa`, non come dimensione. La versione 1.0 eliminava tale
coordinata durante la pulizia dei metadati GRIB; di conseguenza il dataset
finale non conteneva la coordinata `level`.

La versione 1.1.0:

- preserva il livello prima della pulizia delle coordinate ausiliarie;
- lo converte in una dimensione `level` espressa in hPa;
- concatena i singoli file lungo l'asse verticale;
- usa il livello riportato nel nome del file come fallback;
- ordina i livelli dal basso verso l'alto atmosferico, cioè da 1000 a 250 hPa.

Non è necessario riscaricare i GRIB già presenti. È sufficiente sostituire lo
script e rilanciare lo stesso comando.

## Confini regionali nei PNG — v1.2.0

Per aggiungere i confini amministrativi regionali alle mappe singole e al
pannello usare l'opzione `--regional-boundaries`:

```powershell
python pyroconvection_icon2i_meteohub_v1_2.py `
  --run latest `
  --lead-hours 0 3 6 9 12 15 18 21 24 `
  --area 45.2 8.8 41.4 13.5 `
  --regional-boundaries `
  --outdir output_piroconvezione_icon2i
```

Il layer predefinito è Natural Earth `admin_1_states_provinces_lines` a scala
1:10m. Al primo utilizzo Cartopy può scaricarlo nella propria cache.

Per utilizzare un proprio shapefile geografico EPSG:4326:

```powershell
python pyroconvection_icon2i_meteohub_v1_2.py ... `
  --regions-shapefile "C:\dati\confini_regioni_wgs84.shp"
```

Lo shapefile abilita automaticamente la sovrapposizione. I confini vengono
stampati nei PNG singoli e nel pannello, non vengono incorporati nei GeoTIFF.
