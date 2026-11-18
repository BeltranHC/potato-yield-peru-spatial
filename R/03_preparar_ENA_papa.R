# =============================================================================
# 03_preparar_ENA_papa.R — base UA-campaña de papa, 2014–2025
# Entrada : ENA_2014_2024.sav (docente, latin1) + data/ENA_2025/ENA_2025.sav (UTF-8)
# Salida  : data/interim/ena_papa_ua_todas.rds  (todas las UA con papa)
#
# Cambios respecto a R/01_pipeline_ENA_espacial.R (versión anterior):
#   * Clave de UA según el diseño muestral del año:
#       2014–2022: CONGLOMERADO + NSELUA + UA
#       2023–2025: NSEGM + ID_PROD + UA   (CONGLOMERADO/NSELUA vienen vacíos)
#     La versión anterior colapsaba todas las UA de 2023–2024 en una sola.
#   * Factor de expansión: FACTOR (hasta 2022) o FACTOR_PRODUCTOR (desde 2023).
#   * Coordenadas: en 2023–2025 son del SEGMENTO (en 2025, redondeadas a 0,01°,
#     ~1,1 km), no de la parcela -> se marca prec_coord para el análisis puntual.
# =============================================================================
source("R/00_config.R")
N_MAX <- as.numeric(Sys.getenv("N_MAX", "Inf"))   # p. ej. N_MAX=2e5 para pruebas

vars <- c("ANIO","CCDD","CCPP","CCDI","CONGLOMERADO","NSELUA","NSEGM","ID_PROD","UA",
          "ESTRATO","REGION","FACTOR","FACTOR_PRODUCTOR","LONGITUD","LATITUD","ALTITUD",
          "P204_NOM","P204_COD","P209_MES","P209_ANIO","P206_INI_MES","P206_INI_ANIO",
          "P207_FIN_MES","P207_FIN_ANIO",
          "P210_SUP_1","P210_SUP_2","P217_SUP_1","P217_SUP_2",
          "P219_CANT_1","P219_CANT_2","P219_EQUIV_KG",
          "P212","P213","P214","P236","P238","P240","P801","P901","P902","P905","P704","P225",
          "P224A")

# ---- Lectura -----------------------------------------------------------------
leer_ena <- function(f, encoding, n_max = Inf) {
  nm <- names(read_sav(f, n_max = 0, encoding = encoding))
  sel <- nm[toupper(nm) %in% vars]
  t0 <- Sys.time()
  d <- read_sav(f, n_max = n_max, encoding = encoding, col_select = all_of(sel)) |> zap_labels()
  names(d) <- toupper(names(d))
  message("  ", basename(f), ": ", nrow(d), " filas, ", ncol(d), " variables (",
          round(difftime(Sys.time(), t0, units = "mins"), 1), " min)")
  d <- mutate(d, across(everything(), as.character))
  for (v in setdiff(vars, names(d))) d[[v]] <- NA_character_
  vacio <- function(x) is.na(x) | trimws(x) == ""
  d |> mutate(
    c1 = ifelse(vacio(CONGLOMERADO), NSEGM, CONGLOMERADO),
    c2 = ifelse(vacio(NSELUA), ID_PROD, NSELUA),
    FACTOR = ifelse(vacio(FACTOR), FACTOR_PRODUCTOR, FACTOR),
    id_ua = paste(ANIO, key(c1), key(c2), key(UA), sep = "_")) |>
    select(-c1, -c2)
}

# ---- Utilidades ----------------------------------------------------------------
nn  <- function(x) suppressWarnings(as.numeric(x))
dec <- function(d) nn(paste0("0.", ifelse(is.na(d) | d == "", "0", d)))   # "125" -> 0.125
num <- function(ent, d) nn(ent) + dec(d)
si  <- function(x) as.integer(nn(x) == 1)                                   # 1 = Sí, 2 = No
maxna <- function(x) { x <- x[!is.na(x)]; if (length(x)) max(x) else NA }
first_ok <- function(x) { x <- x[!is.na(x) & x != ""]; if (length(x)) x[1] else NA }
fecha <- function(a, m) {
  ok <- !is.na(nn(a)) & !is.na(nn(m)) & nn(m) %in% 1:12
  out <- rep(as.Date(NA), length(a))
  out[ok] <- as.Date(sprintf("%04d-%02d-01", as.integer(nn(a[ok])), as.integer(nn(m[ok]))))
  out
}
ubigeo_de <- function(dd, pp, di) ifelse(is.na(nn(dd)) | is.na(nn(pp)) | is.na(nn(di)), NA,
  sprintf("%02d%02d%02d", as.integer(nn(dd)), as.integer(nn(pp)), as.integer(nn(di))))

# ---- Base UA-campaña de papa ------------------------------------------------------
preparar_ena <- function(ena) {
  geo <- ena |>
    mutate(ubigeo = ubigeo_de(CCDD, CCPP, CCDI)) |>
    group_by(id_ua) |>
    summarise(ANIO = as.integer(first_ok(ANIO)), ubigeo = first_ok(ubigeo),
              REGION = nn(first_ok(REGION)), ESTRATO = first_ok(ESTRATO),
              FACTOR = nn(first_ok(FACTOR)),
              LONGITUD = nn(first_ok(LONGITUD)), LATITUD = nn(first_ok(LATITUD)),
              ALTITUD = nn(first_ok(ALTITUD)), .groups = "drop") |>
    mutate(prec_coord = case_when(is.na(LATITUD) ~ NA_character_,
                                  ANIO >= 2023 ~ "segmento", TRUE ~ "UA"))

  hh <- ena |>
    group_by(id_ua) |>
    summarise(fert = maxna(si(P238)), abono = maxna(si(P236)), plag = maxna(si(P240)),
              asoc = maxna(si(P801)),
              credito = maxna(ifelse(nn(P901) == 2, 0L, si(P902))),   # P902 solo si P901 = 1
              seguro = maxna(si(P905)), asist = maxna(si(P704)),
              horas_cap = maxna(nn(P225)), .groups = "drop")

  papa <- ena |>
    # "^PAPA" a secas también captura PAPAYA y PAPAYUELA: exigir "PAPA" sola o "PAPA <variedad>"
    filter(grepl("^PAPA($|[[:space:]])", toupper(trimws(P204_NOM))), !is.na(nn(P219_EQUIV_KG))) |>
    transmute(id_ua, cultivo = P204_NOM,
      sup_sem   = num(P210_SUP_1, P210_SUP_2),
      sup_cos   = num(P217_SUP_1, P217_SUP_2),
      prod_kg   = num(P219_CANT_1, P219_CANT_2) * nn(P219_EQUIV_KG),  # cantidad x kg/unidad
      riego     = as.integer(nn(P212) %in% 2:7),                     # 1 = lluvia (secano)
      riego_tec = as.integer(nn(P213) %in% 1:4),
      semilla_cert = as.integer(nn(P214) == 1),
      siembra = fecha(P209_ANIO, P209_MES),
      cosecha = fecha(P207_FIN_ANIO, P207_FIN_MES)) |>
    filter(sup_cos > 0, prod_kg > 0)

  papa |>
    group_by(id_ua) |>
    summarise(n_cultivos = n(), prod_kg = sum(prod_kg), sup_cos = sum(sup_cos),
              sup_sem = sum(sup_sem, na.rm = TRUE),
              riego = max(riego), riego_tec = max(riego_tec), semilla_cert = maxna(semilla_cert),
              siembra = suppressWarnings(min(siembra, na.rm = TRUE)),
              cosecha = suppressWarnings(max(cosecha, na.rm = TRUE)), .groups = "drop") |>
    mutate(rend = prod_kg / sup_cos, y = log(rend), ln_sup = log(sup_cos),
           perdida = pmax(0, 1 - sup_cos / pmax(sup_sem, sup_cos))) |>
    filter(rend >= 500, rend <= 60000) |>
    inner_join(geo, by = "id_ua") |>
    left_join(hh, by = "id_ua")
}

# ---- Ejecución -------------------------------------------------------------------
message("ENA 2014–2024 (docente) ...")
papa_ua <- preparar_ena(leer_ena(f_ena_docente, "latin1", N_MAX)); invisible(gc())

if (file.exists(f_ena25)) {
  message("ENA 2025 ...")
  papa_ua <- bind_rows(papa_ua, preparar_ena(leer_ena(f_ena25, "UTF-8", N_MAX)))
} else message("Sin ", f_ena25, ": ejecutar R/02_armar_ENA2025.R")

saveRDS(papa_ua, file.path(dir_interim, "ena_papa_ua_todas.rds"))

resumen <- papa_ua |> group_by(ANIO) |>
  summarise(n_ua = n(), n_distritos = n_distinct(ubigeo), rend_med_t = median(rend) / 1000,
            riego = mean(riego), credito = mean(credito, na.rm = TRUE),
            con_coord = mean(!is.na(LATITUD)), horas_cap_NA = mean(is.na(horas_cap)))
print(resumen, n = Inf, width = Inf)
write_csv(resumen, file.path(dir_tab, "T0_resumen_papa_por_anio.csv"))
