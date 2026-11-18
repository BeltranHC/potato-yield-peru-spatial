# Spatially non-stationary drivers of potato yield in the Peruvian Andes (ENA 2015–2025)

**English summary.** This repository contains the R and Python code for the article
*"Spatially non-stationary drivers of potato yield in the Peruvian Andes: a multiscale and
geostatistical analysis of the National Agricultural Survey (2015–2025)"*. It reproduces:
- the ESDA (Moran's I, LISA with FDR, variograms);
- the spatial econometric models (SAR, SEM, SDM with impacts);
- MGWR with exact inference;
- the cross-scale consistency analysis;
- regression kriging under random, segment-grouped, spatial-block, kNNDM and leave-year-out
  cross-validation, with an external 2025 test.

The raw data are not included. ENA microdata are public at the INEI microdata portal
(https://proyectos.inei.gob.pe/microdatos/), and the covariates (TerraClimate, SRTM, travel
time) are downloaded by `R/04_covariables.R`. Aggregated result tables (`output/tablas`) and
figures (`output/figuras`, `output/figuras_en`) are included. Run the scripts in numerical
order from the RStudio project root; script 10 also requires Python with `pip install mgwr`.
Comments in the code are in Spanish.

---

# artfred — Rendimiento de papa en los Andes peruanos (ENA 2015–2025)

Proyecto RStudio del artículo. Abrir `artfred.Rproj`: el directorio de trabajo es esta carpeta.
Todos los scripts empiezan con `source("R/00_config.R")`.

## Scripts (ejecutar en orden)

| Script | Qué hace | Estado |
|---|---|---|
| `R/00_config.R` | Rutas, paquetes, constantes (CRS, URL INEI) | ✔ |
| `R/01_descargar_ENA2025.R` | Descarga y descomprime los 21 módulos ENA 2025 (1893–1913) | ✔ |
| `R/02_armar_ENA2025.R` | Arma `data/ENA_2025/ENA_2025.sav` con la estructura del docente | ✔ |
| `R/03_preparar_ENA_papa.R` | Base UA-campaña de papa 2014–2025 | ✔ |
| `R/04_covariables.R` | Límites, SRTM 90 m, tiempo de viaje, TerraClimate (ppt, tmin, tmax, déficit) → UA y distritos | ✔ |
| `R/05_ESDA.R` | Tabla 1, Moran global (sensibilidad a W y por año), LISA-FDR, variograma | ✔ |
| `R/06_modelos_distritales.R` | OLS → LM → SLX/SAR/SEM/SDM/SDEM + impactos; GWR y MGWR | ✔ |
| `R/06b_robustez_riego.R` | Riego: interacción con déficit, EF distrito/segmento, tecnificado vs. gravedad | ✔ |
| `R/funciones_prediccion.R` | Datos y funciones comunes de predicción (las usan 07 y 09) | ✔ |
| `R/08_revision_consistencia.R` | Revisión 1 y 5: consistencia entre escalas y robustez del SDM (7 variantes, incluye panel) | ✔ |
| `R/09_revision_prediccion.R` | Revisión 2 y 4: CV por segmento-año y por ubicación, calibración, smearing, bootstrap | ✔ |
| `R/10_mgwr_inferencia.R` + `py/mgwr_inferencia.py` | Revisión 3: MGWR con ENP exacto, IC del ancho de banda y número de condición local (requiere `pip install mgwr`) | ✔ |
| `R/11_figuras_en.R` | Versiones en inglés de las figuras (`output/figuras_en/`) | ✔ |
| `R/12_contraste_oficial.R` | Contraste de la ENA expandida con las cifras oficiales de MIDAGRI | ✔ |
| `R/07_RK_validacion.R` | OLS/RK/RF/RFRK; CV aleatoria, bloques, kNNDM y por año; prueba 2025; mapa con AOA | ✔ |

## Manuscrito

`manuscrito/manuscript_en.qmd` es la versión en inglés para el envío; usa las figuras de
`output/figuras_en` (generadas con `R/11_figuras_en.R`). `manuscrito/manuscrito.qmd` (en español) lee las tablas de `output/tablas` y las figuras
de `output/figuras`. Si se vuelven a correr los scripts, las tablas del artículo se actualizan
solas; las cifras del texto están escritas a mano y hay que revisarlas. Para exportar a Word,
desde RStudio usar el botón *Render*, o en una terminal:

    quarto render manuscrito/manuscrito.qmd --to docx

Pendiente en el texto: nombres de autores, cifras oficiales de MIDAGRI en la introducción,
enlace al repositorio y verificar las referencias en Scopus.

## Datos

* `../ENA_2014_2024/ENA_2014_2024.sav`: base del docente (14,7 GB). No se copia; se
  lee con `encoding = "latin1"`.
* `data/` no se versiona en git (ver `.gitignore`).
  * `data/raw/ENA_2025/`: zip y .sav por módulo + cuestionarios y diccionarios PDF.
  * `data/ENA_2025/ENA_2025.sav`: se lee con `encoding = "UTF-8"`. Los CSV `diag_*`
    documentan la armonización.

## Hallazgos sobre la ENA (para la sección de Datos del artículo)

* **Cambio de diseño muestral en 2023.** La clave de UA pasa de
  `CONGLOMERADO + NSELUA + UA` a `NSEGM + ID_PROD + UA`, y el factor de expansión
  pasa de `FACTOR` a `FACTOR_PRODUCTOR`. 2025 sigue el esquema 2023–2024.
* **Coordenadas.** En 2023–2025 son del segmento, no de la parcela. En 2025 vienen
  redondeadas a 0,01° (~1,1 km), lo que implica un error de posición ≥ 1 km en el
  análisis puntual.
* **Filtro de papa.** Con `^PAPA` también entraban PAPAYA y PAPAYUELA (652 + 23 registros
  solo en 2025), lo que daba falsos clústeres en la selva. Ahora el filtro exige "PAPA" sola
  o "PAPA <variedad>": blanca, amarilla, color, nativa, huayro, amarga.
* **ESDA (script 05).** Moran global del ln-rendimiento: I = 0,43 (kNN-6), 0,37 (kNN-10) y
  0,48 (reina), todos con p < 0,001. Por año, I va de 0,24 a 0,37, siempre significativo.
  LISA con FDR: 44 distritos Alto-Alto y 48 Bajo-Bajo de 778. Variograma de UA: pepita/meseta
  ≈ 0,68 y rango práctico del residuo ≈ 40 km, que es el tamaño de bloque para la validación.
* **Riego (script 06b).** El signo negativo en el distrito es una falacia ecológica. Dentro
  del mismo distrito, el riego aumenta el rendimiento: +3,5 % (EF distrito) y +4,4 % (EF
  segmento×año). El efecto crece con el déficit hídrico (+6 % por cada DE, p < 0,001). El
  riego tecnificado da +10,7 % y el riego por gravedad no tiene efecto significativo (−2,4 %).
  Ponderado por FACTOR: −2,7 % (n.s.).
* **Validación (script 07).** Entrenamiento con 17 211 UA (2023–2024) y prueba externa con
  9 461 UA (2025).
  * La CV aleatoria infla el NSE entre 0,16 y 0,24 frente a la prueba 2025. kNNDM y bloques
    quedan a ±0,03 de esa prueba.
  * Mejor modelo: RFRK con manejo, NSE = 0,24 en 2025. Con solo ambiente (mapeable), RK y
    RFRK llegan a ≈ 0,15.
  * Para el kriging, los residuos se promedian por coordenada única (hay segmentos
    compartidos) y la pepita se trata como error de medición.
  * El 61 % de la rejilla de 5 km queda dentro del AOA.
  * `F7_mapa_RK_2025.png` está desactualizada (el archivo estaba bloqueado al escribir).
    Usar `F7_mapa_RK_2025_AOA.png`.
* **Covariables (script 04).**
  * Clima: TerraClimate (público, ~4 km, 2013–2025), en la ventana siembra–cosecha de cada
    UA. El 83 % de las UA declara una ventana válida. En 2021 no hay fechas, así que se usa
    oct–abr. PISCOp/ERA5 quedan como robustez, porque requieren cuenta.
  * Acceso: tiempo de viaje a ciudades de ≥ 20 mil habitantes (Nelson et al., 2019).
    Reemplaza a P225, que solo existe hasta 2019.
  * Distritos: 1 566 con papa y 865 con ≥ 15 UA. Hay 34 UBIGEO nuevos, reasignados a la
    capital de su provincia.
  * Control de calidad: el 93 % de los puntos 2023–2025 cae en su distrito declarado.
* **Estructura de la base del docente.** Dentro de cada UA, registros de módulos
  distintos comparten fila. Las variables de cabecera de la UA se repiten en todas
  las filas.
