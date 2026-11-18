# =============================================================================
# 00_config.R — rutas, paquetes y utilidades comunes del proyecto artfred
# Se ejecuta al inicio de cada script con: source("R/00_config.R")
# Directorio de trabajo = raíz del proyecto RStudio (artfred/)
# =============================================================================

options(timeout = 3600, scipen = 999, encoding = "UTF-8")

# ---- Paquetes ----------------------------------------------------------------
pk_base <- c("haven", "dplyr", "tidyr", "purrr", "stringr", "readr")
pk_esp  <- c("sf", "terra", "exactextractr", "spdep", "spatialreg",
             "GWmodel", "gstat", "CAST", "blockCV", "ggplot2")

instalar_faltantes <- function(pk) {
  falta <- setdiff(pk, rownames(installed.packages()))
  if (length(falta)) install.packages(falta)
  invisible(falta)
}
instalar_faltantes(pk_base)
suppressPackageStartupMessages(invisible(lapply(pk_base, library, character.only = TRUE)))

# ---- Rutas -------------------------------------------------------------------
# La base del docente (14,7 GB) NO se copia: se lee desde la carpeta padre.
f_ena_docente <- "../ENA_2014_2024/ENA_2014_2024.sav"

dir_data    <- "data"
dir_raw25   <- file.path(dir_data, "raw", "ENA_2025")       # zip y .sav por módulo
dir_zip25   <- file.path(dir_raw25, "zip")
dir_sav25   <- file.path(dir_raw25, "sav")
dir_ena25   <- file.path(dir_data, "ENA_2025")              # salida: ENA_2025.sav
f_ena25     <- file.path(dir_ena25, "ENA_2025.sav")
dir_interim <- file.path(dir_data, "interim")               # .rds intermedios
dir_ext     <- file.path(dir_data, "ext")                   # PISCO, ERA5, SRTM, MODIS, límites
dir_out     <- "output"
dir_tab     <- file.path(dir_out, "tablas")
dir_fig     <- file.path(dir_out, "figuras")

for (d in c(dir_zip25, dir_sav25, dir_ena25, dir_interim, dir_ext, dir_tab, dir_fig))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- ENA 2025 (INEI, código de encuesta 1036) --------------------------------
url_ena25   <- "https://proyectos.inei.gob.pe/iinei/srienaho/descarga/SPSS/1036-Modulo%d.zip"
modulos25   <- 1893:1913
modulos_clave25 <- c(`1893` = "Características de la UA (I)",
                     `1894` = "Características de la UA (II) / ubicación / FACTOR",
                     `1895` = "Producción agrícola y riego",
                     `1896` = "Pérdidas",
                     `1899` = "Insumos",
                     `1907` = "Asistencia técnica / extensión",
                     `1908` = "Asociatividad",
                     `1909` = "Crédito y seguro",
                     `1911` = "Productor")

# Clave de UA: ANIO + CONGLOMERADO + NSELUA + UA ("0001" == "1")
key <- function(x) {
  n <- suppressWarnings(as.integer(as.character(x)))
  ifelse(is.na(n), as.character(x), as.character(n))
}

# Límites distritales limpios (st_make_valid deja GEOMETRYCOLLECTION y vacías)
leer_distritos <- function(f = file.path(dir_ext, "limites", "distritos.gpkg")) {
  d <- sf::st_read(f, quiet = TRUE)
  d <- d[!sf::st_is_empty(d), ]
  suppressWarnings(sf::st_collection_extract(d, "POLYGON")) |>
    group_by(ubigeo, distrito, provincia, departamento) |> summarise(.groups = "drop") |>
    sf::st_cast("MULTIPOLYGON")
}

CRS_GEO <- 4326
CRS_M   <- 32718   # UTM 18S
