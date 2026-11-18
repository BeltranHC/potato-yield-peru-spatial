# =============================================================================
# 02_armar_ENA2025.R — arma data/ENA_2025/ENA_2025.sav con la MISMA estructura
# que la base del docente (ENA_2014_2024.sav): mismas 1 679 columnas, mismo orden,
# mismos tipos y etiquetas, ANIO = 2025 y módulos combinados en filas por UA.
#
# Cómo está armada la base del docente (verificado sobre 2023–2024):
#   * Desde 2023 la clave de UA es ANIO + NSEGM + ID_PROD + UA (antes:
#     CONGLOMERADO + NSELUA + UA) y el factor es FACTOR_PRODUCTOR (antes FACTOR).
#   * No es un apilado puro: dentro de cada UA, registros de módulos distintos
#     comparten fila (2024: 845 mil filas frente a ~2,5 M de registros de módulo).
#   * Las variables de cabecera de la UA (ubicación, estrato, factor...) se
#     repiten en todas las filas.
# Regla usada aquí (reproducible y sin pérdida de información):
#   1. "Cabecera" = columnas que aparecen en >= 2 módulos y son constantes en la UA.
#   2. Los módulos se reparten en "carriles": dos módulos van en el mismo carril
#      solo si no comparten ninguna columna que no sea de cabecera; así ningún valor
#      se sobrescribe. Dentro del carril, el k-ésimo registro de cada módulo en la UA
#      ocupa la fila k. Los carriles se apilan.
#   3. Se completan las columnas del docente ausentes en 2025 con NA.
# Salidas: data/ENA_2025/ENA_2025.sav y diagnósticos CSV en data/ENA_2025/.
# =============================================================================
source("R/00_config.R")

# ---- 0. Armonización de nombres ----------------------------------------------
# Si el INEI renombró una variable en 2025: NOMBRE_DOCENTE = "NOMBRE_2025".
# Revisar data/ENA_2025/diag_sugerencias_renombre.csv (emparejado por etiqueta).
renombrar_2025 <- c()

clave25 <- c("ANIO", "NSEGM", "ID_PROD", "UA")

# ---- 1. Plantilla del docente (solo cabecera; ~2 s) --------------------------
f_tmpl <- file.path(dir_interim, "plantilla_ENA_docente.rds")
if (!file.exists(f_tmpl)) {
  saveRDS(read_sav(f_ena_docente, n_max = 0, encoding = "latin1"), f_tmpl)
}
tmpl <- readRDS(f_tmpl)
TN <- names(tmpl); TU <- toupper(TN)
# Las etiquetas de la plantilla vienen con doble codificación ("AÃ±o"): se reparan.
fix_moji <- function(s) {
  if (is.null(s)) return(s)
  r <- iconv(s, "UTF-8", "latin1"); Encoding(r) <- "UTF-8"
  ifelse(is.na(r) | !validUTF8(r) | !grepl("[ÂÃ]", s), s, r)
}
etq <- function(x) { l <- attr(x, "label", exact = TRUE); if (is.null(l)) NA_character_ else fix_moji(l) }

# ---- 2. Lectura de módulos 2025 (vienen en UTF-8: NO usar encoding latin1) ---
fs <- list.files(dir_sav25, "[.]sav$", recursive = TRUE, full.names = TRUE)
if (!length(fs)) stop("Primero ejecute R/01_descargar_ENA2025.R")
names(fs) <- sub("[.]sav$", "", basename(fs), ignore.case = TRUE)

leer_mod <- function(f) {
  d <- read_sav(f)
  names(d) <- toupper(names(d))
  if (length(renombrar_2025)) {
    ren <- setNames(toupper(renombrar_2025), toupper(names(renombrar_2025)))
    d <- rename(d, any_of(ren))
  }
  d
}
message("Leyendo ", length(fs), " archivos de módulo 2025 ...")
mods <- lapply(fs, leer_mod)

# ---- 3. Diagnóstico de columnas ----------------------------------------------
cols25 <- bind_rows(lapply(names(mods), function(m)
  tibble(archivo = m, variable = names(mods[[m]]),
         etiqueta = vapply(mods[[m]], etq, ""), filas = nrow(mods[[m]]))))
nuevas <- cols25 |> filter(!variable %in% TU)
faltan <- tibble(variable = TN, etiqueta = vapply(tmpl, etq, "")) |>
  filter(!toupper(variable) %in% cols25$variable)

norm <- function(s) tolower(gsub("[^[:alnum:]]", "", iconv(s, "UTF-8", "ASCII//TRANSLIT")))
sug <- inner_join(mutate(faltan, k = norm(etiqueta)),
                  mutate(nuevas, k = norm(etiqueta)) |> select(k, variable_2025 = variable, archivo),
                  by = "k", relationship = "many-to-many") |> filter(k != "")

write_csv(cols25, file.path(dir_ena25, "diag_columnas_2025.csv"))
write_csv(nuevas, file.path(dir_ena25, "diag_nuevas_2025_no_en_docente.csv"))
write_csv(faltan, file.path(dir_ena25, "diag_docente_ausentes_en_2025.csv"))
write_csv(sug,    file.path(dir_ena25, "diag_sugerencias_renombre.csv"))
message("Columnas del docente presentes en 2025: ", length(TN) - nrow(faltan), " de ", length(TN),
        " | nuevas 2025 (se descartan): ", n_distinct(nuevas$variable),
        " | sugerencias de renombre: ", nrow(sug))

# ---- 4. Llevar cada módulo a los tipos/etiquetas de la plantilla -------------
coercion_log <- list()
a_plantilla <- function(x, t, v, m) {
  x0 <- zap_labels(x)
  if (is.character(t)) {
    y <- as.character(x0)
  } else {
    y <- suppressWarnings(as.numeric(x0))
    perd <- sum(is.na(y) & !is.na(x0) & trimws(as.character(x0)) != "")
    if (perd) coercion_log[[length(coercion_log) + 1]] <<- tibble(archivo = m, variable = v, na_por_coercion = perd)
  }
  vl <- attr(x, "labels", exact = TRUE)
  if (is.null(vl)) { vl <- attr(t, "labels", exact = TRUE); if (!is.null(vl)) names(vl) <- fix_moji(names(vl)) }
  lab <- etq(t)
  if (!is.null(vl)) {
    vl <- if (is.character(y)) setNames(as.character(vl), names(vl)) else setNames(as.numeric(vl), names(vl))
    vl <- vl[!duplicated(vl)]
    y <- labelled(y, vl, label = lab)
  } else if (!is.na(lab)) attr(y, "label") <- lab
  y
}
for (m in names(mods)) {
  d <- mods[[m]][, names(mods[[m]]) %in% TU, drop = FALSE]
  for (v in names(d)) d[[v]] <- a_plantilla(d[[v]], tmpl[[match(v, TU)]], v, m)
  if (!all(clave25 %in% names(d))) stop("Módulo sin clave de UA: ", m)
  d$.k <- paste(key(d$NSEGM), key(d$ID_PROD), key(d$UA), sep = "_")
  mods[[m]] <- d
}
coercion_log <- bind_rows(coercion_log)
if (nrow(coercion_log)) {
  write_csv(coercion_log, file.path(dir_ena25, "diag_na_por_coercion.csv"))
  warning(nrow(coercion_log), " variables con valores perdidos al convertir a numérico (ver diag_na_por_coercion.csv)")
}

# ---- 5. Cabecera de UA: columnas compartidas y constantes dentro de la UA ----
n_mod <- table(unlist(lapply(mods, function(d) setdiff(names(d), ".k"))))
compartidas <- setdiff(names(n_mod)[n_mod >= 2], clave25)
es_constante <- function(v) {
  z <- bind_rows(lapply(mods, function(d) if (v %in% names(d))
    tibble(.k = d$.k, val = as.character(zap_labels(d[[v]])))))
  z <- z |> filter(!is.na(val), val != "") |> distinct()
  !anyDuplicated(z$.k)
}
cabecera <- c(clave25, Filter(es_constante, compartidas))
no_cab   <- setdiff(compartidas, cabecera)
message("Cabecera de UA (", length(cabecera), "): ", paste(cabecera, collapse = " "))
if (length(no_cab)) message("Compartidas NO constantes (fuerzan carriles distintos): ",
                            paste(no_cab, collapse = " "))

cab <- bind_rows(lapply(mods, function(d) d[, intersect(c(".k", cabecera), names(d))])) |>
  group_by(.k) |>
  summarise(across(everything(), ~ { z <- .x[!is.na(.x) & as.character(.x) != ""]; if (length(z)) z[1] else .x[1] }),
            .groups = "drop")

# ---- 6. Carriles: módulos sin columnas propias en común comparten filas ------
# Los módulos de nivel parcela se alinean por número de parcela (P105): los datos de
# una misma parcela de módulos distintos quedan en la misma fila.
alinear <- "P105"
propias <- lapply(mods, function(d) setdiff(names(d), c(".k", cabecera, alinear)))
orden   <- names(sort(sapply(mods, nrow), decreasing = TRUE))
carril  <- list()
for (m in orden) {
  ok <- which(vapply(carril, function(cm) !length(intersect(propias[[m]], unlist(propias[cm]))), TRUE))
  if (length(ok)) carril[[ok[1]]] <- c(carril[[ok[1]]], m) else carril[[length(carril) + 1]] <- m
}
asig <- bind_rows(lapply(seq_along(carril), function(i) tibble(n_carril = i, archivo = carril[[i]])))
print(asig, n = Inf)

armar_carril <- function(cm) {
  partes <- lapply(cm, function(m) {
    d <- mods[[m]][, c(".k", propias[[m]]), drop = FALSE]
    p <- if (alinear %in% names(mods[[m]])) as.character(zap_labels(mods[[m]][[alinear]])) else NA
    d$.p <- ifelse(is.na(p), "", p)
    d$.s <- ave(seq_len(nrow(d)), d$.k, d$.p, FUN = seq_along)
    d
  })
  Reduce(function(a, b) full_join(a, b, by = c(".k", ".p", ".s")), partes)
}
message("Combinando ", length(carril), " carriles ...")
ena25 <- bind_rows(lapply(seq_along(carril), function(i) mutate(armar_carril(carril[[i]]), .carril = i)))
ena25[[alinear]] <- a_plantilla(ifelse(ena25$.p == "", NA, ena25$.p), tmpl[[match(alinear, TU)]], alinear, "alinear")
ena25 <- left_join(ena25, cab, by = ".k")

# Control de integridad: para cada variable, los valores no vacíos de salida deben
# ser la suma de los no vacíos de todos los módulos que la contienen.
nz <- function(x) sum(!is.na(x) & as.character(zap_labels(x)) != "")
v_ctrl <- setdiff(unique(unlist(lapply(mods, names))), c(".k", cabecera, alinear))
# P105 es la llave de alineación: se verifica que no se pierda ningún par UA–parcela
par_in  <- unique(unlist(lapply(mods, function(d) if (alinear %in% names(d))
  paste(d$.k, zap_labels(d[[alinear]]))[!is.na(d[[alinear]])])))
par_out <- unique(paste(ena25$.k, zap_labels(ena25[[alinear]]))[!is.na(ena25[[alinear]])])
if (!setequal(par_in, par_out)) stop("Pares UA-parcela distintos entre entrada y salida")
ctrl <- tibble(variable = v_ctrl,
  modulos = vapply(v_ctrl, function(z) paste(names(mods)[vapply(mods, function(d) z %in% names(d), TRUE)], collapse = "+"), ""),
  no_vacios_entrada = vapply(v_ctrl, function(z) sum(vapply(mods, function(d) if (z %in% names(d)) nz(d[[z]]) else 0, 0)), 0),
  no_vacios_salida  = vapply(v_ctrl, function(z) nz(ena25[[z]]), 0))
write_csv(ctrl, file.path(dir_ena25, "diag_control_modulos.csv"))
if (any(ctrl$no_vacios_entrada != ctrl$no_vacios_salida)) {
  print(filter(ctrl, no_vacios_entrada != no_vacios_salida)); stop("Se perdieron valores al combinar")
}
message("Control de integridad OK (", nrow(ctrl), " variables verificadas)")
rm(mods); invisible(gc())

# ---- 7. Estructura final = plantilla del docente -----------------------------
ena25$ANIO <- a_plantilla(rep(2025, nrow(ena25)), tmpl$ANIO, "ANIO", "constante")
ena25 <- ena25 |> arrange(CCDD, CCPP, CCDI, NSEGM, ID_PROD, UA, .carril, .p, .s)
n <- nrow(ena25)
falt_cols <- setdiff(TU, names(ena25))
for (v in falt_cols) {
  t <- tmpl[[match(v, TU)]]
  ena25[[v]] <- a_plantilla(if (is.character(t)) rep(NA_character_, n) else rep(NA_real_, n), t, v, "vacía")
}
ena25 <- ena25[, TU]
names(ena25) <- TN                                       # mayúsculas/minúsculas del docente
stopifnot(identical(names(ena25), names(tmpl)),
          all(vapply(seq_along(TN), function(i) is.character(ena25[[i]]) == is.character(tmpl[[i]]), TRUE)))

# ---- 8. Resumen ----------------------------------------------------------------
message("Filas ENA_2025: ", n, " | UA: ", n_distinct(paste(ena25$NSEGM, ena25$ID_PROD, ena25$UA)),
        " | columnas: ", ncol(ena25))

# ---- 9. Guardar ---------------------------------------------------------------
message("Escribiendo ", f_ena25, " ...")
write_sav(ena25, f_ena25, compress = "byte")
message("Listo: ", f_ena25, " (", round(file.size(f_ena25) / 2^30, 2), " GB). ",
        "Leer con haven::read_sav(..., encoding = 'UTF-8').")
