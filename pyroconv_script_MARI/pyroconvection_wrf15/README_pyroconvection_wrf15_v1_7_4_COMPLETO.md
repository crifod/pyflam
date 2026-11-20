# Piroconvezione con WRF 1.5 km ItaliaMeteo–CIMA–LaMMA

Documentazione complessiva dello script:

```text
pyroconvection_wrf15_v1_7_4.py
```

Lo script elabora i campi del modello operativo **WRF a 1,5 km** distribuiti tramite MeteoHub e produce una classificazione della **predisposizione atmosferica alla piroconvezione**.

La versione 1.7.4 conserva integralmente la ladder adattiva, le soglie e la correzione ABL della v1.7.3. La novità riguarda la selezione temporale: `--start` e `--end` non sono più obbligatori.

Con il download automatico è quindi possibile elaborare direttamente l'ultima corsa disponibile: 

```powershell
python pyroconvection_wrf15_v1_7_4.py --download-wrf --wrf-run latest
```

Quando `--start` e `--end` sono entrambi omessi, lo script usa tutte le scadenze valide disponibili nella corsa selezionata, mantenendo soltanto le ore UTC indicate da `--hours`. Se si desidera limitare l'intervallo, i due parametri devono essere specificati insieme.

La correzione ABL introdotta nella v1.7.3 interviene esclusivamente quando si verificano contemporaneamente queste condizioni:

1. il primo livello originale disponibile sopra 200 m supera già `RiB = 0,33`, quindi il crossing non è delimitato;
2. `PBLH/HPBL` non è finito o non è positivo;
3. il primo livello sopra soglia si trova a quota non superiore a 600 m AGL.

In questo caso lo script **non inventa un'altezza ABL**: conserva soltanto il limite superiore fisicamente dimostrabile e assegna classe 0 `Surface plume`. Se il limite superiore è maggiore di 600 m, la cella resta NoData. Tutti gli altri step orari seguono esattamente la logica della v1.7.2.

La classificazione combina:

- altezza dello strato limite atmosferico, ABL;
- altezza del livello di condensazione, LCL;
- stabilità media nell'ABL, misurata sui livelli originali oppure rappresentata tramite un proxy;
- gradiente di temperatura potenziale sopra l'ABL;
- quota del massimo wind shear, solo quando realmente risolta;
- umidità relativa in prossimità del top dell'ABL;
- modalità, qualità e motivazione della classificazione.

> **Avvertenza operativa**  
> Le classi descrivono una predisposizione atmosferica sotto l'ipotesi di un incendio sufficientemente intenso. Non rappresentano una probabilità calibrata di PyroCb e devono essere validate su casi locali. La classe va interpretata insieme a `classification_mode`, `classification_quality` e `classification_reason`.

---

## 1. Dataset WRF

Lo script è predisposto per il dataset operativo WRF ItaliaMeteo–CIMA–LaMMA disponibile su MeteoHub.

Caratteristiche considerate:

- previsione deterministica;
- risoluzione orizzontale circa 1,5 km;
- campi di superficie e livelli isobarici;
- corse tipiche 00 e 12 UTC;
- orizzonte previsionale fino a circa +72 ore;
- file NetCDF, GRIB, GRIB2, GRB o GRB2.

Il programma può:

1. scaricare automaticamente i file open data;
2. utilizzare file già scaricati;
3. distinguere automaticamente i file di superficie da quelli isobarici;
4. ricostruire la dimensione verticale `level` quando i livelli sono distribuiti in file GRIB separati;
5. gestire griglie rettilinee e griglie WRF curvilinee con coordinate bidimensionali;
6. applicare `--bbox` e `--stride` ai singoli file prima del merge;
7. valutare la risoluzione verticale sui livelli originali;
8. selezionare automaticamente la ladder `full`, `resolved_noshear` o `proxy_noshear`;
9. riconoscere una ABL certamente poco profonda tramite un limite superiore RiB, senza attribuirle una quota puntuale;
10. distinguere classe 0 e NoData;
11. mascherare il mare mediante LSM WRF o tramite fallback geografico;
12. disegnare per impostazione predefinita i confini provinciali;
13. usare automaticamente l’intero intervallo valido della corsa quando `--start` e `--end` sono omessi.

---

## 2. Classi, modalità, qualità e motivazione

### 2.1 Classi atmosferiche: `pyro_class`

| Valore | Classe | Significato sintetico |
|---:|---|---|
| 0 | Surface plume | Plume di superficie; ABL troppo poco sviluppata o condizioni insufficienti per una colonna convettiva classificabile |
| 1 | Convective plume | Colonna convettiva senza PyroCu persistente |
| 2 | Overshooting PyroCu | Potenziale superamento dell'ABL con PyroCu breve |
| 3 | Resilient PyroCu | PyroCu persistente, ma crescita verticale limitata |
| 4 | PyroCb | Predisposizione a deep PyroCu/PyroCb |
| 255 | NoData | Profilo non classificabile, dati insufficienti o punto escluso |

Il valore 0 rappresenta una **classe fisica valida**. Il valore 255 rappresenta invece una cella non classificabile.

### 2.2 Modalità adattiva: `classification_mode`

| Valore | Modalità | Significato |
|---:|---|---|
| 0 | `invalid` | Profilo non classificabile |
| 1 | `proxy_noshear` | Stabilità ABL non risolta direttamente; uso del proxy della profondità di mixing e ladder senza shear obbligatorio |
| 2 | `resolved_noshear` | Stabilità ABL misurata sui livelli originali, ma ladder completa non disponibile |
| 3 | `full` | Stabilità, cap, RH al top e shear risolti sui livelli originali |
| 4 | `shallow_abl_upper_bound` | Nessuna altezza ABL puntuale: il primo crossing RiB è certamente entro 600 m e basta per assegnare Surface plume |

La modalità viene scelta per ogni cella e per ogni tempo valido. La modalità 4 non è una ladder termodinamica completa: è un'inferenza conservativa a qualità bassa.

### 2.3 Qualità: `classification_quality`

| Valore | Qualità | Interpretazione |
|---:|---|---|
| 0 | Non valida | NoData/non classificabile |
| 1 | Bassa | Proxy, fallback PBLH o campionamento verticale limitato |
| 2 | Media | Diagnostica principale risolta, ma ladder completa non disponibile |
| 3 | Alta | Ladder completa e diagnostiche verticali risolte |

La qualità viene ridotta quando:

- l'ABL proviene da `PBLH/HPBL` anziché da un crossing Richardson risolto;
- il cap sopra l'ABL non è risolto;
- la RH al top ABL non è adeguatamente racchiusa dai livelli originali;
- la classe 0 deriva soltanto da un limite superiore RiB, senza una quota ABL puntuale (`classification_mode = 4`).

### 2.4 Motivazione: `classification_reason`

Il nuovo campo rende auditabile la ladder e permette di capire perché una cella è bianca, verde, gialla, arancione o rossa.

| Codice | Significato | Classe associata normalmente |
|---:|---|---:|
| 0 | `invalid_or_nodata` | 255 |
| 1 | `shallow_abl_surface_plume` | 0 |
| 2 | `stable_or_poorly_mixed_convective_plume` | 1 |
| 3 | `lcl_too_high_convective_plume` | 1 |
| 4 | `overshooting_pyrocu` | 2 |
| 5 | `insufficient_deep_mixing_convective_plume` | 1 |
| 6 | `dry_abl_top_convective_plume` | 1 |
| 7 | `rh_top_unresolved_resilient_pyrocu` | 3 |
| 8 | `cap_unresolved_resilient_pyrocu` | 3 |
| 9 | `cap_not_weak_resilient_pyrocu` | 3 |
| 10 | `shear_not_close_resilient_pyrocu` | 3 |
| 11 | `deep_pyrocu_noshear` | 4 |
| 12 | `deep_pyrocu_full` | 4 |
| 13 | `shallow_abl_upper_bound_surface_plume` | 0 |

Questo campo è particolarmente utile per analizzare le ore notturne e distinguere:

- classe 0 dovuta ad ABL inferiore a 600 m;
- classe 1 dovuta a strato stabile;
- classe 1 dovuta a LCL troppo elevato;
- classe 1 dovuta a top ABL troppo secco;
- classe 0 ottenuta da un limite superiore RiB certo, senza trasformare indiscriminatamente i NoData in Surface plume.

---

## 3. Requisiti software

Versione Python consigliata:

```text
Python 3.10 o superiore
```

Installazione con `pip`:

```powershell
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install numpy pandas xarray netCDF4 requests matplotlib rasterio cfgrib eccodes shapely
```

Per produrre mappe con coste, confini nazionali e regionali:

```powershell
pip install cartopy
```

Conda:

```powershell
conda create -n pyroconvection python=3.11
conda activate pyroconvection
conda install -c conda-forge numpy pandas xarray netcdf4 requests matplotlib rasterio cfgrib eccodes shapely cartopy
```

Verifica:

```powershell
python pyroconvection_wrf15_v1_7_4.py --help
```

---

## 4. Variabili richieste

### 4.1 Campi di superficie

- temperatura a 2 m;
- temperatura di rugiada a 2 m oppure umidità relativa a 2 m;
- vento U e V a 10 m;
- PBLH/HPBL;
- orografia o geopotenziale superficiale;
- pressione al suolo oppure pressione ridotta al livello del mare;
- land-sea mask, se disponibile.

Alias principali:

```text
2t, t2m, T2
2d, d2m, TD2
10u, u10, U10
10v, v10, V10
hpbl, blh, PBLH
orog, HGT, terrain_height, z_sfc
sp, PSFC, surface_pressure
prmsl, msl, MSLP
lsm, LANDMASK, landmask, land_sea_mask
```

Se la pressione al suolo non è presente, viene stimata da PRMSL, orografia e T2m.

### 4.2 Campi ai livelli isobarici

Sono richiesti:

- geopotenziale o altezza geopotenziale;
- temperatura;
- umidità specifica oppure umidità relativa;
- vento U e V.

Shortname MeteoHub tipici:

```text
gh t q u v
```

### 4.3 Livelli predefiniti

```text
1000 975 950 925 900 850 800 750 700 650
600 550 500 450 400 350 300 hPa
```

Modifica:

```powershell
--pressure-levels 1000 950 925 900 850 800 700 600 500 400 300
```

---

## 5. Download automatico da MeteoHub

### 5.1 Esempio operativo Toscana: ultima corsa, intervallo automatico

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --download-wrf `
  --wrf-run latest `
  --download-dir dati_wrf `
  --hours 00 03 06 09 12 15 18 21 `
  --bbox 44.6 9.7 42.2 12.6 `
  --stride 2 `
  --class-only `
  --regional-boundaries `
  --outdir output_piroconvezione_wrf
```

In questo caso lo script individua la corsa più recente e usa automaticamente tutte le sue scadenze disponibili. Il filtro `--hours` resta attivo: con i valori dell’esempio vengono elaborate le scadenze ogni tre ore.

Per una corsa specifica:

```powershell
--wrf-run 2026071400
```

Ordine del bounding box:

```text
NORTH WEST SOUTH EAST
```

### 5.2 Variabili scaricate per impostazione predefinita

```text
2t 2d 10u 10v hpbl orog prmsl lsm gh t q u v
```

Il campo `lsm` viene richiesto, ma se non è disponibile lo script può derivare la maschera terrestre dai confini provinciali.

### 5.3 Struttura dei download

```text
dati_wrf\<RUN>\<VARIABILE>\<FILE>
```

Viene creato anche:

```text
download_manifest.txt
```

Il downloader supporta:

- file temporanei `.part`;
- ripresa del download;
- mancato riscaricamento dei file completi;
- `--overwrite-downloads`;
- segnalazione delle variabili mancanti.

---

## 6. Uso di file locali

### 6.1 File misti

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --wrf-files dati_wrf\2026071400\*\*.grib `
  --start 2026-07-14 `
  --end 2026-07-16 `
  --hours 00 03 06 09 12 15 18 21 `
  --bbox 44.6 9.7 42.2 12.6 `
  --outdir output_piroconvezione_wrf
```

### 6.2 File locali con intervallo automatico

Anche con file già presenti è possibile omettere le date. Vengono usate tutte le scadenze contenute nei file, filtrate da `--hours`:

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --wrf-files dati_wrf\2026071400\*\*.grib `
  --hours 00 03 06 09 12 15 18 21 `
  --bbox 44.6 9.7 42.2 12.6 `
  --outdir output_piroconvezione_wrf
```

### 6.3 File surface e pressure separati

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --surface-files dati\surface\*.grib `
  --pressure-files dati\pressure\*.grib `
  --start 2026-07-14 `
  --end 2026-07-16 `
  --bbox 44.6 9.7 42.2 12.6 `
  --outdir output_piroconvezione_wrf
```

### 6.4 NetCDF multidimensionale

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --wrf-files dati\WRF_run00.nc `
  --start 2026-07-14 `
  --end 2026-07-16 `
  --outdir output_piroconvezione_wrf
```

---

## 7. Tempo valido e ore UTC

Il filtro usa il **tempo valido della previsione**, non soltanto l'ora di inizializzazione.

Sono gestite:

- `time`: run di riferimento nei GRIB;
- `step`: scadenza;
- `valid_time`: tempo valido;
- `forecast_reference_time`: run conservato nell'output normalizzato.

Sono disponibili due modalità.

### 7.1 Intervallo automatico

Omettendo entrambi i parametri, lo script utilizza il primo e l’ultimo tempo valido realmente presenti nei file o nella corsa scaricata:

```powershell
--download-wrf --wrf-run latest --hours 00 03 06 09 12 15 18 21
```

Nel log viene riportato, per esempio:

```text
Intervallo temporale automatico dalla corsa: 2026-07-15 00:00 -> 2026-07-18 00:00 UTC (...)
```

### 7.2 Intervallo esplicito

Per limitare le date, `--start` e `--end` devono essere indicati insieme:

```powershell
--start 2026-07-14 --end 2026-07-16 --hours 00 03 06 09 12 15 18 21
```

Non è ammesso specificare soltanto uno dei due estremi.

Le ore sono UTC. Durante l'ora legale italiana:

```text
00Z = 02:00 locali
03Z = 05:00 locali
06Z = 08:00 locali
18Z = 20:00 locali
21Z = 23:00 locali
```

Questa conversione è importante nell'interpretazione del ciclo diurno dell'ABL.

---

## 8. Ottimizzazione di velocità e memoria

La catena applica:

1. apertura del singolo file;
2. normalizzazione di variabili, coordinate, livello e tempo;
3. selezione anticipata delle scadenze;
4. applicazione del bounding box al singolo file;
5. applicazione dello stride;
6. merge/concat dei dataset già ridotti;
7. selezione dei livelli;
8. allineamento delle griglie;
9. elaborazione di un timestamp alla volta.

### 8.1 Bounding box pre-merge

```powershell
--bbox 44.6 9.7 42.2 12.6
```

### 8.2 Stride

```powershell
--stride N
```

| Stride | Celle conservate indicativamente | Passo indicativo da 1,5 km |
|---:|---:|---:|
| 1 | 100% | 1,5 km |
| 2 | 25% | 3 km |
| 3 | 11% | 4,5 km |
| 4 | 6% | 6 km |

### 8.3 Indici cfgrib persistenti

```powershell
--cfgrib-index-dir C:\lavoro\cache_cfgrib
```

### 8.4 Modalità `--class-only`

Con questa opzione vengono sempre conservati:

```text
pyro_class
classification_mode
classification_quality
classification_reason
land_sea_mask
```

Le diagnostiche verticali estese non vengono allocate.

---

## 9. Griglia, mare e confini amministrativi

### 9.1 Griglie supportate

Rettangolari:

```text
latitude(latitude)
longitude(longitude)
```

Curvilinee WRF:

```text
latitude(y, x)
longitude(y, x)
```

oppure:

```text
XLAT(south_north, west_east)
XLONG(south_north, west_east)
```

### 9.2 Maschera terra/mare

La sequenza è:

1. uso di `lsm`, `LANDMASK`, `landmask` o `land_sea_mask` se presenti;
2. normalizzazione automatica da 0–1 o 0–100;
3. fallback mediante poligoni delle province italiane se LSM è assente;
4. assegnazione di 255 alle celle marine.

Per non mascherare il mare:

```powershell
--keep-sea
```

### 9.3 Confini provinciali

Sono disegnati per impostazione predefinita nei PNG e nel pannello.

Il file viene scaricato nella cache:

```text
.admin_boundaries\limits_IT_provinces.geojson
```

File locale alternativo:

```powershell
--provinces-file C:\dati\confini\province.geojson
```

Disattivazione:

```powershell
--no-province-boundaries
```

Confini regionali aggiuntivi:

```powershell
--regional-boundaries
```

---

## 10. Metodo fisico e classificazione adattiva

### 10.1 Preparazione del profilo

Per ogni cella e tempo:

1. ricostruzione delle quote AGL;
2. eliminazione dei livelli sotto il terreno e sopra 8.000 m;
3. richiesta di almeno quattro livelli validi;
4. calcolo di theta e theta virtuale;
5. stima ABL;
6. calcolo LCL con Bolton;
7. valutazione di stabilità, cap, RH e shear;
8. selezione della ladder;
9. assegnazione di classe, modalità, qualità e motivazione.

L'interpolazione a 50 m è usata per alcune diagnostiche continue, ma **non viene considerata nuova informazione verticale**.

### 10.2 ABL da Bulk Richardson

La soglia è:

```text
RiB >= 0,33 sopra 200 m AGL
```

#### Controllo di risoluzione introdotto nella v1.7.2

Il crossing viene accettato solo quando esistono:

- un livello originale inferiore con `RiB < 0,33`;
- un livello originale superiore con `RiB >= 0,33`.

Il punto di crossing viene interpolato soltanto fra questi due livelli.

Se il **primo livello disponibile sopra 200 m** supera già la soglia, il crossing non è considerato risolto. Il superamento potrebbe infatti essere avvenuto molto più in basso, soprattutto:

- di notte;
- nelle zone collinari e montane;
- quando 1000 hPa è sotto il terreno;
- quando il primo livello disponibile è 975 o 950 hPa.

In questo caso lo script usa `PBLH/HPBL` del WRF.

#### Fallback conservativo aggiunto nella v1.7.3

Al primo step della corsa può accadere che `PBLH/HPBL` sia `NaN`, zero o non ancora inizializzato. La v1.7.2 lasciava correttamente la cella come NoData, producendo però ampie aree grigie anche quando il profilo conteneva un'informazione sufficiente a escludere le classi 1–4.

La v1.7.3 distingue quattro casi:

| Caso | Azione |
|---|---|
| Crossing RiB delimitato | Usa l'altezza interpolata fra i due livelli originali |
| Crossing non delimitato e PBLH valida | Usa PBLH, senza altre modifiche alla ladder |
| Crossing non delimitato, PBLH non valida, primo superamento ≤600 m | Non stima ABL; assegna classe 0 usando soltanto il limite superiore |
| Crossing non delimitato, PBLH non valida, primo superamento >600 m | Mantiene NoData perché non è dimostrabile che ABL sia sotto la soglia operativa |

La conclusione del terzo caso è rigorosa rispetto alla definizione adottata: se il primo livello già sopra soglia è, per esempio, a 480 m, il primo crossing `RiB = 0,33` deve trovarsi a quota minore o uguale a 480 m e quindi certamente sotto 600 m. Lo script non salva `480 m` come altezza ABL: salva `abl_height = NaN` e `abl_upper_bound = 480 m`.

Questa diramazione è isolata e non cambia i risultati degli step nei quali PBLH è valida o il crossing Richardson è delimitato.

Codifica:

| Campo | Valore | Significato |
|---|---:|---|
| `abl_source` | 0 | ABL non disponibile |
| `abl_source` | 1 | crossing Bulk Richardson risolto |
| `abl_source` | 2 | PBLH/HPBL fallback |
| `abl_source` | 3 | solo limite superiore RiB; nessuna altezza ABL puntuale |
| `abl_upper_bound` | quota in m | primo livello sopra soglia usato come limite superiore, non come stima ABL |
| `rib_crossing_resolved` | 0 | crossing non risolto o non usato |
| `rib_crossing_resolved` | 1 | crossing delimitato da livelli originali |

### 10.3 Filtro ABL a 600 m

La v1.7.2 e la v1.7.3 applicano:

```text
ABL < 600 m → classe 0 Surface plume
ABL >= 600 m → prosecuzione della ladder
```

Questa modifica riduce la sovraclassificazione notturna della classe 1. Nella v1.7.3 la stessa regola può essere applicata anche senza conoscere l'altezza esatta, ma soltanto quando `abl_upper_bound <= 600 m` rende la conclusione certa.

La soglia è definita in:

```python
min_classifiable_abl_m = 600.0
```

### 10.4 Stabilità ABL

Il gradiente viene misurato sui livelli originali fra circa il 20% e l'80% dell'ABL.

È considerato risolto con almeno:

```text
3 livelli originali
250 m di estensione verticale
```

Soglie:

```text
dtheta/dz <= 3,0e-3 K m-1 → compatibile almeno con classe 2
dtheta/dz <= 1,0e-3 K m-1 → compatibile con classi persistenti/deep
```

Se non è risolto viene usato il proxy:

```text
ABL >= 455 m
```

Il proxy non supera il filtro principale di 600 m: una ABL di 500 m resta quindi classe 0.

### 10.5 Cap sopra l'ABL

Layer:

```text
ABL + 200 m → ABL + 1200 m
```

Risolto con almeno due livelli originali e 350 m di estensione.

Soglie:

```text
gamma-theta <= 4,2e-3 K m-1 → cap debole
gamma-theta >= 4,8e-3 K m-1 → cap forte
```

Un cap non risolto non abilita la classe 4.

### 10.6 RH al top ABL

Media in:

```text
ABL - 150 m → ABL + 150 m
```

Requisito per classi persistenti:

```text
RH >= 80%
```

Il top deve essere racchiuso da livelli originali con gap non superiore a 900 m.

### 10.7 Massimo shear

Viene calcolato esclusivamente fra livelli originali, indicativamente fra 200 e 6.000 m AGL.

È risolto se:

- sono presenti almeno cinque livelli;
- il massimo non è sul primo o ultimo intervallo;
- i gap locali non superano 1.000 m.

In modalità `full`, la classe 4 richiede:

```text
min(|zshear - ABL|, |zshear - LCL|) / ABL <= 0,30
```

### 10.8 Sequenza della ladder

1. profilo insufficiente → 255;
2. ABL <600 m → classe 0, reason 1;
3. mixing non compatibile con PyroCu → classe 1, reason 2;
4. LCL/ABL >1,60 → classe 1, reason 3;
5. 1,00 ≤ LCL/ABL ≤1,60 → classe 2, reason 4;
6. mixing insufficiente per sviluppo persistente → classe 1, reason 5;
7. RH top <80% → classe 1, reason 6;
8. RH non risolta → classe 3, reason 7;
9. cap non risolto → classe 3, reason 8;
10. cap non debole → classe 3, reason 9;
11. shear lontano in modalità full → classe 3, reason 10;
12. deep ladder noshear → classe 4, reason 11;
13. deep ladder full → classe 4, reason 12.

Lo score 0–100 resta separato e non modifica la classe.

---

## 11. Output

```text
output_piroconvezione_wrf\
├── pyroconvection_wrf15.nc
├── pyroconvection_panel.png
├── classi_piroconvezione.csv
├── geotiff\
│   ├── pyro_class_YYYYMMDDTHHMMZ.tif
│   └── ...
└── png\
    ├── pyro_class_YYYYMMDDTHHMMZ.png
    └── ...
```

### 11.1 Campi sempre presenti

| Variabile | Descrizione |
|---|---|
| `pyro_class` | Classe 0–4; 255 NoData |
| `classification_mode` | Ladder utilizzata |
| `classification_quality` | Qualità 0–3 |
| `classification_reason` | Motivo determinante della classe |
| `land_sea_mask` | Maschera/frazione terra |

### 11.2 Diagnostiche del NetCDF completo

| Variabile | Descrizione | Unità/codifica |
|---|---|---|
| `pyro_score` | Score diagnostico | 0–100 |
| `abl_height` | ABL | m |
| `lcl_height` | LCL | m |
| `lcl_abl_ratio` | Rapporto LCL/ABL | 1 |
| `maximum_shear_height` | Quota massimo shear | m |
| `shear_distance_ratio` | Distanza normalizzata shear | 1 |
| `theta_gradient_abl` | Gradiente theta ABL | K m⁻¹ |
| `gamma_theta_free_troposphere` | Gradiente theta sopra ABL | K m⁻¹ |
| `relative_humidity_abl_top` | RH top ABL | % |
| `abl_source` | Origine o stato diagnostico ABL | 0/1/2/3 |
| `abl_upper_bound` | Limite superiore del crossing RiB non risolto; non è una stima ABL | m |
| `rib_crossing_resolved` | Crossing RiB verticalmente risolto | 0/1 |
| `stability_method` | Proxy o gradiente originale | 0/1/2 |
| `n_original_levels_valid` | Livelli validi | numero |
| `n_original_levels_abl` | Livelli nel layer ABL | numero |
| `n_original_levels_cap` | Livelli nel cap | numero |
| `n_original_levels_shear` | Livelli nel layer shear | numero |
| `cap_resolved` | Cap risolto | 0/1 |
| `rh_top_resolved` | RH top risolta | 0/1 |
| `shear_resolved` | Shear risolto | 0/1 |

### 11.3 GeoTIFF

- `uint8`;
- classi 0–4;
- NoData 255;
- DEFLATE;
- EPSG:4326;
- trasformazione affine per griglie rettilinee;
- GCP per griglie curvilinee.

---

## 12. Opzioni principali della riga di comando

| Opzione | Descrizione |
|---|---|
| `--download-wrf` | Download automatico |
| `--wrf-run latest` | Ultima corsa |
| `--wrf-run YYYYMMDDHH` | Corsa specifica |
| `--download-dir DIR` | Directory download |
| `--download-variables ...` | Variabili MeteoHub |
| `--wrf-files ...` | File misti/locali |
| `--surface-files ...` | File di superficie |
| `--pressure-files ...` | File isobarici |
| `--pressure-levels ...` | Livelli utilizzati |
| `--cfgrib-index-dir DIR` | Cache indici cfgrib |
| `--bbox N W S E` | Ritaglio pre-merge |
| `--stride N` | Sottocampionamento pre-merge |
| `--start` / `--end` | Intervallo valido facoltativo; devono essere usati insieme. Se omessi, viene usato l’intero intervallo della corsa |
| `--hours ...` | Ore UTC |
| `--keep-sea` | Non maschera il mare |
| `--class-only` | Salva i cinque campi essenziali |
| `--no-province-boundaries` | Non disegna le province |
| `--provinces-file FILE` | Confini provinciali locali |
| `--regional-boundaries` | Aggiunge confini regionali |
| `--regions-shapefile FILE` | Confini regionali locali |
| `--no-individual-png` | Disabilita PNG singoli |
| `--no-panel-png` | Disabilita pannello |
| `--panel-cols N` | Colonne pannello |
| `--no-cartopy` | Grafici lat/lon semplici |
| `--verbose` | Log dettagliato |

---

## 13. Esempi operativi

### 13.1 Produzione completa

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --download-wrf `
  --wrf-run latest `
  --hours 00 03 06 09 12 15 18 21 `
  --bbox 44.6 9.7 42.2 12.6 `
  --stride 1 `
  --regional-boundaries `
  --outdir output_piroconvezione_wrf
```

### 13.2 Produzione rapida

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --download-wrf `
  --wrf-run latest `
  --hours 00 06 12 18 `
  --bbox 44.6 9.7 42.2 12.6 `
  --stride 2 `
  --class-only `
  --outdir output_piroconvezione_wrf_rapido
```

### 13.3 Test molto rapido

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --wrf-files dati\*.grib `
  --start 2026-07-14 `
  --end 2026-07-14 `
  --hours 00 06 12 18 `
  --bbox 44.6 9.7 42.2 12.6 `
  --stride 4 `
  --class-only `
  --outdir output_test_stride4
```

---

## 14. Messaggi di log attesi

```text
File WRF separati automaticamente: ... superficie, ... isobarici
Filtro livelli isobarici: apro ... file e ne escludo ...
Ritaglio pre-merge: griglia ...x... -> ...x..., stride=...
Intervallo temporale automatico dalla corsa: ... -> ... UTC (... scadenze)
Livelli isobarici disponibili: ... hPa
Dominio: ... tempi, ... righe, ... colonne, ... livelli
Maschera terra/mare: uso il campo ... presente nei dati WRF
```

oppure:

```text
Campo LSM/LANDMASK assente: land-sea mask derivata geometricamente da ...
```

Al termine:

```text
Classificazione adattiva: validi=...; proxy-noshear=..., resolved-noshear=..., full=..., shallow-upper-bound=...
Motivi principali: ABL calcolata < 600 m=...; ABL certamente < 600 m da limite RiB=...; mixed layer stabile/poco rimescolato=...
```

Quest'ultima riga consente un controllo immediato della prevalenza notturna della classe 0 e della classe 1.

---

## 15. Risoluzione dei problemi

### 15.1 Prima mappa grigia o molte celle 255 al passo +0 h

Controllare nel NetCDF:

```text
pyro_class
classification_mode
classification_reason
abl_height
abl_upper_bound
abl_source
rib_crossing_resolved
```

Interpretazione della v1.7.3:

- `classification_reason = 1`: ABL puntuale calcolata o ottenuta da PBLH e inferiore a 600 m;
- `classification_reason = 13`, `classification_mode = 4`, `abl_source = 3`: crossing non delimitato, PBLH non valida e limite superiore entro 600 m; cella correttamente assegnata a Surface plume;
- `abl_height = NaN` con `abl_upper_bound <= 600 m` e motivo 13 è intenzionale: la quota esatta non è nota;
- `pyro_class = 255` con `abl_upper_bound > 600 m`: NoData mantenuto correttamente, perché il profilo non dimostra che l'ABL sia sotto 600 m;
- `pyro_class = 255` e `abl_upper_bound = NaN`: verificare livelli, temperatura, umidità, vento, pressione e geopotenziale.

La v1.7.3 non trasforma genericamente il grigio in bianco. Interviene soltanto sul sottoinsieme di celle per cui il primo crossing è matematicamente vincolato entro 600 m. Pertanto una parte residua di grigio può essere corretta e deve restare visibile.

### 15.1.1 Molto verde nelle ore notturne

Controllare:

```text
classification_reason
classification_quality
abl_source
```

- motivo 2: strato stabile o poco rimescolato;
- motivo 3: LCL troppo alto;
- motivo 6: top ABL troppo secco;
- `abl_source = 2`: uso di PBLH WRF.

Dalla v1.7.2, ABL inferiori a 600 m restano classe 0. La v1.7.3 estende questa protezione ai soli casi con limite superiore certo. Un residuo di verde è plausibile quando PBLH o un crossing delimitato indicano ABL almeno pari a 600 m.

### 15.2 Molti fallback PBLH

Può essere normale di notte o sui rilievi. Verificare:

```text
abl_source = 2
rib_crossing_resolved = 0
```

Il fallback evita di assegnare artificialmente come ABL la quota del primo livello isobarico sopra 200 m.

### 15.3 Molte celle 255

Verificare:

- almeno quattro livelli validi;
- pressione al suolo coerente;
- disponibilità di T, q/RH, U, V e geopotenziale;
- maschera terra/mare;
- date e ore richieste.

### 15.4 Apertura lenta o memoria elevata

Usare:

1. `--bbox`;
2. `--stride 2` o `--stride 3`;
3. `--class-only`;
4. meno ore;
5. livelli strettamente necessari;
6. indici cfgrib persistenti.

### 15.5 Confini provinciali non disponibili

Specificare un file locale:

```powershell
--provinces-file C:\dati\confini\province.geojson
```

### 15.6 GeoTIFF curvilineo

Regolarizzazione:

```powershell
gdalwarp -t_srs EPSG:4326 -r near input.tif output.tif
```

---

## 16. Differenze rispetto alla procedura ICON-2I

La versione WRF:

- usa fino a 17 livelli isobarici anziché 5–6;
- tenta la ladder completa solo quando la risoluzione originale lo consente;
- passa automaticamente a `noshear` quando necessario;
- conserva il proxy della profondità di mixing;
- applica la soglia ABL di 600 m;
- accetta l'ABL Richardson solo con crossing verticalmente delimitato;
- espone motivazione e qualità per ogni cella;
- non applica il gate dei combustibili e rappresenta quindi il solo potenziale atmosferico.

---

## 17. Evoluzione delle versioni

### Versione 1.0–1.4

- conversione iniziale da ERA5 a WRF;
- download MeteoHub;
- normalizzazione del tempo valido;
- griglie curvilinee;
- ricostruzione della dimensione `level`;
- separazione surface/pressure.

### Versione 1.5

- selezione anticipata delle ore e dei livelli;
- indici cfgrib persistenti;
- elaborazione di un timestamp alla volta;
- `float32` e `--class-only`.

### Versione 1.6

- `--bbox` applicato ai file prima del merge;
- introduzione di `--stride`;
- riuso degli indexer spaziali.

### Versione 1.7

- ladder adattiva `full`, `resolved_noshear`, `proxy_noshear`;
- controlli sui livelli originali;
- classe 0 separata da 255;
- qualità e modalità;
- shear non più inventato sulla griglia a 50 m;
- eliminazione della promozione basata sullo score.

### Versione 1.7.1

- maschera mare da LSM WRF;
- fallback geografico con poligoni provinciali;
- confini provinciali nei PNG e nel pannello;
- supporto `--provinces-file` e cache amministrativa;
- `lsm` aggiunto al download automatico.

### Versione 1.7.2

- soglia classificabile ABL portata da 300 a 600 m;
- `ABL <600 m` forzata a classe 0 Surface plume;
- nuovo crossing Bulk Richardson che richiede un vero bracket verticale;
- fallback PBLH quando il primo livello sopra 200 m supera già `RiB = 0,33`;
- nuovo campo `rib_crossing_resolved`;
- nuovo campo sempre presente `classification_reason`;
- codifica esplicita dei motivi 0–12;
- log riepilogativo delle celle con ABL poco profonda e mixed layer stabile;
- `--class-only` esteso alla motivazione della classe.

### Versione 1.7.3

- correzione mirata del primo step WRF quando PBLH è nulla o non inizializzata;
- conservazione dell'altezza del primo superamento RiB come **limite superiore**, non come stima ABL;
- assegnazione conservativa della classe 0 soltanto se tale limite è ≤600 m;
- mantenimento del NoData quando il limite è >600 m o non disponibile;
- nuova modalità `classification_mode = 4`, `shallow_abl_upper_bound`;
- nuovo motivo `classification_reason = 13`;
- nuovo codice `abl_source = 3`;
- nuova diagnostica `abl_upper_bound`;
- risultati invariati rispetto alla v1.7.2 quando PBLH è valida o il crossing RiB è delimitato;
- log separato delle celle classificate tramite limite superiore.

---

## 18. Uso scientifico e operativo

Il prodotto deve essere interpretato insieme a:

- intensità e potenza potenziale dell'incendio;
- disponibilità e continuità del combustibile;
- vento al suolo e in quota;
- profilo di umidità;
- instabilità atmosferica;
- osservazioni radar, satellitari e fulminazioni;
- evoluzione reale della colonna convettiva.

La classe 4 non implica che un PyroCb si formerà. Indica che il profilo atmosferico è favorevole allo sviluppo verticale profondo qualora il forcing termico dell'incendio sia sufficiente.

Per confrontare elaborazioni diverse devono rimanere coerenti:

- versione dello script;
- soglie `Thresholds`;
- livelli isobarici;
- bbox e stride;
- ore richieste;
- disponibilità del PBLH e delle variabili verticali.

---

## 19. Riferimento metodologico

Castellnou et al. (2022), quadro concettuale sulla transizione fra plume convettivo, PyroCu e PyroCb.

La parametrizzazione resta prototipale e richiede:

- confronto con radiosondaggi;
- verifica su casi storici positivi e negativi;
- analisi separata giorno/notte;
- calibrazione della soglia ABL di 600 m;
- confronto tra ABL da Richardson e PBLH WRF;
- verifica delle classi 4 ottenute in modalità `noshear`.

---

## 19.1 Novità della versione 1.7.4

- `--start` e `--end` sono facoltativi;
- se entrambi sono omessi, viene selezionato automaticamente l’intero intervallo disponibile;
- `--hours` continua a filtrare le sole ore UTC desiderate anche in modalità automatica;
- se viene indicato un solo estremo temporale, lo script termina con un messaggio esplicito;
- il log riporta l’intervallo temporale effettivamente selezionato;
- rimane pienamente compatibile l’uso tradizionale con entrambe le date.

---

## 20. Comando sintetico consigliato

```powershell
python pyroconvection_wrf15_v1_7_4.py `
  --download-wrf `
  --wrf-run latest `
  --hours 00 03 06 09 12 15 18 21 `
  --bbox 44.6 9.7 42.2 12.6 `
  --stride 2 `
  --class-only `
  --regional-boundaries `
  --cfgrib-index-dir C:\lavoro\cache_cfgrib `
  --outdir output_piroconvezione_wrf
```

Con `--class-only` restano disponibili:

```text
pyro_class
classification_mode
classification_quality
classification_reason
land_sea_mask
```

Per analizzare nel dettaglio le ore notturne, eseguire senza `--class-only` e controllare anche:

```text
abl_height
abl_source
rib_crossing_resolved
stability_method
relative_humidity_abl_top
```
