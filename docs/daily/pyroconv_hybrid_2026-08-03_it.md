---
title: "Previsione del tipo di piroconvezione in Toscana -- livelli modello ICON-EU + filtro ICON-2I 2.2 km -- VALIDA 2026-08-03"
subtitle: "Trioraria, 24 h. Corsa 2026-08-03 00Z. Metodo secondo Castellnou et al. (2022), JGR-Atmos."
geometry: a4paper, landscape, margin=1.2cm
fontsize: 9pt
lang: it
header-includes: |
  \usepackage{float}
  \floatplacement{figure}{H}
---

## Come si legge questo prodotto

Ogni mappa risponde a una domanda. Vanno lette nell'ordine -- ma si veda l'avvertenza sotto la
domanda 3: le domande **non** sono annidate, e una cella puo superare la 3 fallendo la 2.

| # | Domanda | Mappa | Il criterio |
|:--|:--|:--|:--|
| 1 | Qual è la piroconvezione più intensa che questa **colonna** potrebbe sostenere? | POTENZIALE | la scala diagnostica, senza alcuna condizione sul fuoco |
| 2 | C'è abbastanza fuoco per generare un **pennacchio**? | FILTRATO PER COMBUSTIBILE | intensità di Byram >= 10 MW/m sui combustibili .lcp |
| 3 | C'è abbastanza fuoco per generare un **pyroCb** in questa specifica colonna? | MARGINE PFT | potenza totale >= PFT propria di quella colonna |

La domanda 1 è una proprietà della sola atmosfera: nessuna informazione sui combustibili o
sul suolo vi entra, nemmeno per stabilire dove sia definita. Le domande 2 e 3 condizionano
entrambe al fuoco e discendono entrambe dallo stesso campo di intensità di Byram, ma **non**
sono lo stesso criterio e non concordano: la (2) confronta la potenza *per metro di fronte*
con una costante fissa di 10 MW/m, mentre la (3) la converte in potenza *totale* attraverso una
lunghezza di fronte di testa assunta e la confronta con una soglia calcolata per quella colonna.
Le celle superano regolarmente la (2) e falliscono la (3) di ordini di grandezza. La domanda 3
è quella operativamente decisiva, ed è la più difficile da soddisfare.

## Il metodo e i suoi limiti (da leggere prima di usare le classi)

La colonna è classificata a partire da un **profilo verticale reale**, non da
un'atmosfera standard: quote per cella dalle quote dei livelli modello ICON-EU (HHL), la profondità di rimescolamento dal **numero di Richardson
bulk** (primo Rib >= 0,33 sopra i 200 m dal suolo, riferito allo stato a 2 m / 10 m), l'LCL
dalla formula esatta di **Bolton (1980)**, e l'umidità alla sommità dell'ABL dalla RH del
modello. Il gamma-theta del tappo è preso fra ABL+200 m e ABL+1200 m.

La **stabilità dello strato rimescolato** è qui una misura autentica, non un proxy. I livelli
modello nativi ICON-EU collocano ~10 livelli dentro lo strato rimescolato (in questa corsa: 15
in mediana sulla terraferma fra le 12 e le 15Z, quando lo strato è maturo -- i conteggi notturni
e serali sono assai più bassi, ma in quelle ore non si classifica nulla). Sul suolo classificato l'interpolazione è supportata nel **99%** delle colonne al culmine della giornata (cioè quella quota contiene almeno 3 livelli dentro lo strato; le restanti restituiscono un gradiente nan e abbandonano la mappa delle classi, ed è anche ciò che assottiglia i pannelli serali). Il dtheta/dz
dello strato rimescolato è dunque una vera interpolazione ai minimi quadrati su quei livelli.
Validata sui radiosondaggi IGRA (GLA 12Z, intero periodo di registrazione), l'interpolazione sui
livelli modello riduce il bias del gradiente a ~-0,4e-4 K/m (da ~-4,7e-4 sui 5 livelli di
pressione ICON-2I) e il bias dell'ABL da Rib a ~-70 m (da ~-700 m). È questa la ragione per cui
il prodotto ibrido esiste: i dati aperti ICON-2I a 2,2 km non risolvono lo strato rimescolato,
ICON-EU sì.

Il prezzo è la risoluzione orizzontale: l'atmosfera è a 6,5 km (riportata sulla griglia a
2,2 km), mentre il **filtro combustibile conserva i campi superficiali ICON-2I a 2,2 km**, dove
il dettaglio orografico conta. Il gradiente dello strato rimescolato è ora misurato, ma la
soglia di 1,1e-3 K/m ricade ancora dentro la barra d'errore validata (+/- ~2,6e-4): i *conteggi*
delle classi vanno perciò trattati come indicativi, non esatti.

Dopo il tramonto lo strato su cui il gradiente viene interpolato è lo **strato residuo** --
l'aria quasi neutra che lo strato convettivo in decadimento lascia dietro di sé -- anziché il
sottile strato stabile notturno che una particella riferita alla superficie troverebbe. Di
giorno i due coincidono. È questo che consente alle ore serali di portare una previsione; è una
diagnostica più recente del resto della scala, e quelle ore vanno pesate di conseguenza.

**Scala effettivamente usata in questa corsa: `shear`.** È girato il metodo completo a 5 diagnostiche: i livelli modello hanno risolto una quota di massimo shear, sicché la classe 4 ha richiesto in aggiunta che quel massimo si collocasse entro 0,30 ABL dall'ABL/LCL -- una condizione *necessaria* che le scale ridotte omettono. La classe 4 porta qui dunque la propria clausola di shear; resta comunque un allerta a ispezionare la colonna, non una probabilità calibrata.

Le classi esprimono la predisposizione atmosferica **dato un incendio di potenza
sufficiente**; nel pannello filtrato quella potenza è calcolata, non presupposta. Sono soglie
derivate dalla letteratura, non validate localmente.

## 1. Potenziale -- che cosa potrebbe sostenere questa colonna? (limite superiore atmosferico)

![Domanda 1 -- POTENZIALE. Classe di piroconvezione dalla sola atmosfera, nessuna condizione sul fuoco. Limite superiore, non un'attesa.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_potential_2026-08-03_it.png){width=100%}

La scala diagnostica con la condizione sul fuoco completamente disattivata: questa è dunque la
risposta della sola atmosfera. È la piroconvezione più intensa che la colonna sosterrebbe **dato
un incendio in grado di sfruttarla**. Quell'incendio è presupposto, non calcolato: entra
attraverso le soglie calibrate della scala, non attraverso una mappa di combustibili. Questa
mappa serve a vedere dove sia l'atmosfera il vincolo stringente, mai come previsione. Si noti
che il classificatore gira su qualunque colonna utilizzabile, anche sul mare: nella figura il
mare è mascherato, ma i raster `pyroconv_potential_*.tif` sottostanti no, quindi ogni statistica
zonale ricavata direttamente da essi va limitata alla terraferma.

## 2. Filtrato per combustibile -- c'è abbastanza fuoco per un pennacchio?

![Domanda 2 -- FILTRATO PER COMBUSTIBILE. La stessa scala diagnostica, con ogni cella sotto i 10 MW/m di intensità di Byram sui combustibili .lcp toscani forzata alla classe 0. È il prodotto atteso.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_gated_2026-08-03_it.png){width=100%}

La stessa scala diagnostica, sulle stesse colonne, con ogni cella forzata alla classe 0 dove i
combustibili .lcp toscani e la meteorologia superficiale prevista non sostengono **10 MW/m** di
intensità del fronte secondo Byram (Tedim et al. 2018). Il filtro è binario e unidirezionale:
sopra la soglia i combustibili non hanno più alcuna influenza, sotto la soglia la cella riporta
un pennacchio superficiale qualunque cosa dica l'atmosfera. Questa mappa non può quindi mai
superare la mappa del potenziale -- ogni differenza fra le due è una cella riportata a 0 -- ed è
il prodotto operativo atteso.

## 3. Margine PFT -- c'è abbastanza fuoco per un pyroCb in *questa* colonna?

![Domanda 3 -- MARGINE PFT. Potenza dell'incendio rapportata alla PyroCb Firepower Threshold propria di ciascuna colonna, scala logaritmica imperniata sul criterio. I cerchi soddisfano il criterio e superano il filtro combustibile; le croci lo soddisfano solo perché la PFT è degenerata.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_pft_margin_2026-08-03_it.png){width=100%}

Potenza **totale** dell'incendio per cella, divisa per la **PyroCb Firepower Threshold** propria
di quella colonna (Tory & Kepert 2021, eq. 31: `PFT = 0.3 z_fc^2 U_ML dtheta_fc`, in GW). Da 1
in su l'incendio è abbastanza potente da forzare un pyroCb *in quella colonna*; sotto 1 non lo
è, per quanto favorevole sia l'atmosfera. Il colore si imperna esattamente su 1 e la scala è
logaritmica, perché il campo si estende su quattro decadi.

L'intensità di Byram è una potenza per metro di fronte e la PFT è una potenza totale: per
collegarle serve una lunghezza attiva del fronte di testa. Anziché assumerne una, la si ricava
dal comportamento del fuoco pubblicato dal Wildfire Data Portal -- la mediana sui suoi 20
incendi che riportano sia il rapporto di superficie bruciata sia la velocità di avanzamento,
700 m -- con una frazione convettiva pari a 0.7 del calore rilasciato che entra nel
pennacchio (Tory & Kepert app. D). Sono entrambe assunzioni, e il margine scala linearmente con
ciascuna: un fronte di testa di 5 km, che Tory & Kepert indicano per un caso estremo,
moltiplicherebbe ogni margine qui riportato per circa sette. **Il margine va letto come un
ordine di grandezza, non come un numero.**

È questa la ragione per cui la mappa conserva valore nonostante l'avvertenza: sulla maggior
parte del dominio bruciabile il margine non è prossimo a 1, ma due o tre decadi al di sotto,
sicché la conclusione «qui nessun pyroCb» è robusta rispetto all'assunzione sul fronte di testa
in un modo in cui una cella marginale non lo sarebbe. Dove invece il margine si avvicina a 1,
l'assunzione diventa portante e la cella merita un giudizio, non la mappa.

### Le domande 2 e 3 non sono annidate -- si leggano i simboli, non il colore

Una cella può soddisfare il criterio PFT **fallendo** il filtro dei 10 MW/m, e in questa corsa
è quanto accade per la maggior parte di esse. Il meccanismo sta al denominatore: la PFT è
`0.3 z_fc^2 U_ML dtheta_fc`, sicché in una colonna bassa, debolmente tappata e con vento
leggero essa collassa -- a fronte di una mediana di dominio prossima a 90 GW, le celle con
margine >= 1 hanno PFT di pochi GW. Un incendio da 2-8 MW/m «supera» allora una soglia che un
incendio di quelle dimensioni non può fisicamente dirsi aver battuto, perché a quell'intensità
non c'è alcun pyroCu da approfondire.

La mappa distingue i due casi, e solo il primo è una previsione:

* **cerchio nero** -- margine >= 1 *e* filtro combustibile superato. Un autentico candidato pyroCb.
* **croce grigia** -- margine >= 1 ma filtro non superato. La PFT è degenerata, non è stata
  battuta. Vanno letti come un artefatto della forma della soglia a piccoli `z_fc`, non come
  un segnale.

La colonna **anche filtro superato** della tabella che segue è dunque il conteggio operativo.
Prendere invece il conteggio grezzo del margine sovrastimerebbe di parecchie volte l'area
candidata.

| Ora | % griglia con potenza | % margine >= 1 | % margine >= 0,1 | margine max | celle >= 1 | **anche filtro superato** |
|:--|--:|--:|--:|--:|--:|--:|
| 00Z | 0.0 | 0.0 | 0.0 | nan | 0 | **0** |
| 03Z | 0.0 | 0.0 | 0.0 | nan | 0 | **0** |
| 06Z | 0.0 | 0.0 | 0.0 | 0.02 | 0 | **0** |
| 09Z | 25.4 | 0.1 | 2.3 | 1.49 | 2 | **0** |
| 12Z | 25.4 | 0.8 | 7.9 | 27.46 | 30 | **16** |
| 15Z | 25.6 | 0.3 | 6.1 | 1.58 | 9 | **7** |
| 18Z | 6.4 | 0.0 | 0.5 | 0.25 | 0 | **0** |
| 21Z | 0.1 | 0.0 | 0.0 | 0.01 | 0 | **0** |

*% griglia con potenza* è la quota di dominio che ha insieme una colonna classificabile e
combustibile bruciabile; le due colonne del margine sono percentuali **di quella quota**, non
della Toscana. Una percentuale sull'intero dominio sarebbe dominata dalle celle prive di
combustibile, che non potranno mai contribuire, e si muoverebbe per ragioni estranee alla
previsione. *celle >= 1* è il conteggio grezzo dietro la prima percentuale, riportato perché a
questi tassi una percentuale si arrotonda a un valore che sembra nulla quando in realtà si
tratta di una manciata di luoghi precisi. ***anche filtro superato*** applica il criterio
congiunto della sezione precedente ed è il conteggio su cui agire; lo scarto fra esso e
*celle >= 1* è la popolazione a PFT degenerata.

## 4. Disaccoppiamento pirogeno secco -- DIAGNOSTICO (nessuna etichetta di classe)

![Diagnostico di supporto -- DISACCOPPIAMENTO PIROGENO SECCO. fireABL/ABL per un incendio intenso di riferimento. Continuo, nessuna etichetta di classe.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_decoupling_2026-08-03_it.png){width=100%}

Il **rapporto di disaccoppiamento** fireABL / ABL indica quanto in alto un incendio intenso di
riferimento (un eccesso di temperatura di 10 K nel pennacchio, l'estremo intenso delle
misure in-plume GRAF, 0,1-13,1 K) accrescerebbe il proprio strato limite per *solo calore
sensibile*, diviso per l'ABL ambientale. È la controparte **secca** delle mappe umide delle
classi qui sopra (Castellnou et al. 2022; Castellnou Ribau et al. 2024): valori ben superiori a
1 segnalano colonne profonde, calde e secche, in cui un incendio può sfondare e disaccoppiarsi
dalla superficie *anche dove la scala umida assegna punteggi bassi* -- che è esattamente la
situazione che la mappa filtrata sottostima. Come la mappa del potenziale, presuppone un incendio
ovunque -- qui con un flusso di riferimento fisso anziché calcolato -- sicché è un limite
superiore, non una previsione; non si applica alcuna distinzione LCL secco/umido (lo scarto di
+1 km presente in letteratura non è supportato dalle etichette prototipo GRAF). fireABL da
`pyflam.atmosphere.fire_induced_abl_grid`, la particella intersecata con la reale colonna
theta(z); le celle che il classificatore rifiuta (nessuna colonna utilizzabile, o ABL sotto i
150 m) restano bianche anziché essere dipinte con un quoziente nudo.

## 5. Altezza di cima prevista -- l'unico campo validato contro osservazioni

![CIMA DEL PENNACCHIO PREVISTA. La scala dei costi risolta per l'altezza invece che per la potenza, per un incendio dichiarato di 5 GW. Continua, nessuna etichetta di classe.](/Users/cristianofoderi/-softEST/firelab-flammap6_install_0828_2025/pyflam/docs/daily/pyroconv_tuscany_hybrid_plumetop_2026-08-03_it.png){width=100%}

Le tre domande precedenti chiedono *quanto fuoco* costi una data quota. Risolvendo la
stessa disequazione nell'altro verso -- per la quota massima che un dato incendio puo
permettersi -- si risponde alla domanda che un previsore ha davvero, e si ottiene un campo
continuo invece di una classe.

E anche l'unico prodotto qui **verificato contro osservazioni da un capo all'altro**.
Confrontare la cima di un pennacchio con la quota-bersaglio di un gradino mette quel bersaglio
su entrambi i lati del test e misura per lo piu la propria circolarita; risolvere per l'altezza
mette un numero previsto contro un numero misurato e nient'altro. Valutato cosi su 1231
pennacchi digitalizzati dalle immagini stereo MISR su 7 regioni e 11 biomi, senza nulla di
adattato: Spearman **+0,46** rispetto alla cima osservata, **+125 m** di bias mediano, 576 m di
errore assoluto medio e 12,8 % meglio del predire una costante su un errore adimensionale --
davanti alla sola potenza (+0,37) e alla sola profondita dello strato limite (+0,34), e con le
mediane regionali ordinate a rho +0,61.

**La potenza va letta come scenario.** La mappa e disegnata per una potenza *totale* dichiarata
di 5 GW, esplicitata qui perche una potenza totale non si ricava da un'intensita di Byram
senza assumere una lunghezza del fronte di testa, ed e proprio quell'assunzione che questo
filone di lavoro esiste per evitare. Raddoppiare l'incendio di riferimento non raddoppia
l'altezza -- il costo cresce con il cubo della salita -- ma il campo si sposta, e non va letto
come stima per cella di cio che brucera davvero.

**Due limiti da dichiarare.** Il campione di validazione e quello di MISR, quindi e fissato
intorno alle 10:30 solari locali: lo strato limite sta ancora crescendo e gli incendi sono piu
piccoli di quanto saranno nel pomeriggio che questo prodotto prevede. E il campo viene prodotto
solo dove la griglia verticale risolve lo strato su cui si legge la cappa, che in questa corsa
significa la sorgente a livelli modello; un ripiego su livelli di pressione lo lascia vuoto
invece di interpolarlo.

A differenza dei raster di classe, `diag_plume_top_*.tif` e `diag_form_*.tif` sono scritti
**gia mascherati** con il campo `valid` del classificatore: una statistica zonale presa
direttamente su di essi non puo includere colonne che il classificatore ha rifiutato. Gli altri
`diag_*.tif` restano grezzi, come sono sempre stati.

## Scala delle classi (da minima a massima attività piroconvettiva)

| Livello | Colore | Classe | Significato |
|:--:|:--|:--|:--|
| 0 | bianco | Pennacchio superficiale | Pennacchio di fumo galleggiante; nessuno sviluppo nuvoloso significativo. |
| 1 | verde | Pennacchio convettivo | Il pennacchio penetra uno strato rimescolato stabile; condensazione possibile, nessun pyroCu. |
| 2 | giallo | PyroCu overshooting | Pirocumulo di breve durata; base della nube sopra l'altezza di rimescolamento (LCL/ABL > 1). |
| 3 | arancione | PyroCu persistente | Pirocumulo persistente in colonna instabile (LCL/ABL < 1). |
| 4 | rosso scuro | PyroCu profondo / pyroCb | Piroconvezione profonda / pirocumulonembo; il tappo superiore debole lascia approfondire il pennacchio. |

## Soglie di classificazione (scala diagnostica in uso: `shear`)

| Diagnostica | Soglia | Effetto |
|:--|:--|:--|
| dtheta/dz strato risc. (minimi quadrati, dentro lo strato) | > 1,1e-3 K/m (stabile) | Solo pennacchio convettivo -- nessun pyroCu |
| dtheta/dz strato risc. (misurato sui livelli modello) | <= 1,1e-3 K/m | Colonna capace di pyroCu (bias ~-0,4e-4 K/m vs radiosondaggi) |
| Rapporto LCL / ABL | 1,0 -- 1,60 | PyroCu overshooting (breve) |
| Rapporto LCL / ABL | < 1,0 | PyroCu persistente |
| Rapporto LCL / ABL | <= 1,10 (+ condizioni sotto) | Ammissibile per pyroCu profondo / pyroCb |
| Gamma-theta del tappo (ABL+200 m -> ABL+1200 m) | <= 4,2e-3 K/m (tappo debole) | Consente l'approfondimento a pyroCb |
| Gamma-theta del tappo | >= 4,8e-3 K/m (tappo forte) | Inibisce l'approfondimento (al più persistente) |
| RH alla sommità dell'ABL (media, ABL +/- 150 m) | >= 60% | Richiesta per le classi 3 e 4 (era 80%; rilassata per condizioni di fire weather secche) |
| Distanza del massimo di shear / ABL | <= 0,30 | Richiesta per la classe 4 -- **risolta e applicata in tutta questa corsa** |
| Intensità del fronte (filtro combustibile, pannello filtrato) | >= 10 MW/m | Potenza minima per qualsiasi pyroCu (Tedim et al. 2018) |
| Profondità dell'ABL | < 150 m | Non classificabile (profondità implausibile). Nessun filtro a 600 m: non aveva base in Castellnou et al. (2022) e scartava 5 degli 8 incendi di campagna -- rimosso il 26/07/2026 |
| Livelli di pressione utilizzabili | < 4 | Cella non classificata |

## Casi di riferimento (Castellnou et al. 2022, Tabella 1)

Gli eventi etichettati nell'articolo, che ancorano la scala pubblicata a 3 diagnostiche
(`pyflam.pyroconvection_type(ladder="castellnou")`, tuttora predefinita nella libreria e
sottoposta a test di regressione su queste righe). Le scale su profilo usate per questa mappa
sono più severe in cima: richiedono in aggiunta una sommità dell'ABL umida e un LCL vicino ad
essa.

| Caso | Tipo osservato | LCL/ABL | dtheta/dz strato risc. | gamma-theta (700-500) |
|:--|:--|:--:|:--|:--:|
| T21 | Pennacchio convettivo | -- | stabile | -- |
| SCQ32 | PyroCu overshooting | > 1 | neutro/instabile | -- |
| M11 | PyroCu persistente | < 1 | instabile | 4,2e-3 (tappo persistente) |
| SCQ41 | PyroCu, non pyroCb | < 1 | instabile | 5,1e-3 (tappo forte) |
| SCQ51 | PyroCu profondo / pyroCb | < 1 | instabile | 3,9e-3 (tappo debole) |

Forzanti: atmosfera dai livelli modello nativi ICON-EU a 6,5 km (dati aperti DWD, CC-BY), i ~24 livelli più bassi; ~15 qui dentro lo strato rimescolato (terraferma, 12-15Z), riportati sulla griglia ICON-2I a 2,2 km. Campi superficiali e filtro combustibile da ICON-2I 2,2 km (MISTRAL / AgenziaItaliaMeteo). Classificatore: pyflam.pyroconvection_type (ABL da Richardson bulk, LCL di Bolton, dtheta/dz misurato nello strato rimescolato, gamma-theta del tappo, RH alla sommità dell'ABL; nessuna CAPE superficiale).
Filtro: Rothermel + chioma Cruz-2005 sui combustibili .lcp con umidità e vento previsti.
Raster diagnostici orari accompagnano le classi: ABL, parcel_ml, residual_ml, LCL, LCL/ABL,
dtheta/dz dello strato rimescolato, gamma-theta del tappo, RH alla sommità, fireABL, decoupling
e (solo ibrido) delta_theta, firecape, penetration, pft_gw (PyroCb Firepower Threshold, GW),
z_fc, delta_theta_fc, u_ml, firepower_gw, pft_margin.
Confini provinciali: derivati ISTAT (openpolis geojson-italy).
Generato da tests/pyroconv_daily.py con pyflam `894a4b9+local-changes`, 2026-08-03 10:21 UTC -- congelato con la fisica di quel
commit, poiché DWD ritira la corsa dopo ~24 h.
