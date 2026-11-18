# =============================================================================
# funciones_prediccion.R — datos y utilidades comunes de la etapa predictiva
# Las usan 07_RK_validacion.R y 09_revision_prediccion.R (no ejecutar solo).
# Requiere: source("R/00_config.R") y los paquetes sf, gstat y ranger cargados.
# =============================================================================
K <- 10; NTREE <- 500; NMAX <- 50

vgi  <- readRDS(file.path(dir_interim, "variograma_ua.rds"))
vm0  <- vgi$vm_res

A <- c("alt_pt", "pend_pt", "ln_tviaje", "ppt_pt", "tmin_pt", "def_pt", "meses_helada_pt")
B <- c(A, "riego", "riego_tec", "semilla_cert", "fert", "abono", "plag", "asist",
       "credito", "cred_na", "asoc", "ln_sup")

# UA 2023–2025 con covariables puntuales, en UTM (m) y con X/Y en km
preparar_puntos <- function(papa) {
  pt <- papa |>
    filter(ANIO >= 2023, !is.na(alt_pt), !is.na(LATITUD)) |>
    mutate(ln_tviaje = log1p(tviaje_pt), cred_na = as.integer(is.na(credito)),
           across(c(credito, semilla_cert, fert, abono, plag, asist, asoc), ~ ifelse(is.na(.x), 0, .x)),
           segmento = sub("^([0-9]+_[^_]+)_.*$", "\\1", id_ua)) |>      # ANIO_NSEGM
    filter(if_all(all_of(c("y", B)), ~ is.finite(.x))) |>
    st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = CRS_GEO) |> st_transform(CRS_M)
  xy <- st_coordinates(pt) / 1000; pt$X <- xy[, 1]; pt$Y <- xy[, 2]
  pt$ubicacion <- paste(round(pt$X, 3), round(pt$Y, 3))
  pt
}

a_km <- function(d) st_as_sf(st_drop_geometry(d), coords = c("X", "Y"), remove = FALSE)
a_err <- function(vm) {                    # Nug -> Err (error de medición)
  m <- as.character(vm$model); m[m == "Nug"] <- "Err"
  vm$model <- factor(m, levels = levels(vgm()$short)); vm$range[m == "Err"] <- 0
  vm
}
krigear <- function(trk, tek, z) {
  vg <- variogram(as.formula(paste(z, "~ 1")), trk, cutoff = 150, width = 5)
  vm <- tryCatch(suppressWarnings(fit.variogram(vg, vm0)), error = function(e) vm0)
  if (any(vm$psill < 0) || any(vm$range[vm$model != "Nug"] <= 0) || isTRUE(attr(vm, "singular"))) vm <- vm0
  # Varias UA comparten la coordenada del segmento -> matriz singular. Se krigea la media
  # de los residuos por ubicación única (el variograma sí usa todas las UA).
  ag <- st_drop_geometry(trk) |> group_by(X, Y) |> summarise(z = mean(.data[[z]]), .groups = "drop") |>
    st_as_sf(coords = c("X", "Y"), remove = FALSE)
  out <- krige(z ~ 1, ag, tek, model = a_err(vm), nmax = NMAX, debug.level = 0)
  list(pred = out$var1.pred, var = out$var1.var, vm = vm)
}
# con_var = TRUE agrega la varianza de predicción de una UA nueva:
#   var_kriging (señal suavizada) + pepita (error de medición, Err)
ajustar_predecir <- function(train, test, vars, con_var = FALSE) {
  f  <- reformulate(vars, "y")
  dtr <- st_drop_geometry(train); dte <- st_drop_geometry(test)
  m  <- lm(f, dtr); p_ols <- predict(m, dte)
  rf <- ranger(f, data = dtr[c("y", vars)], num.trees = NTREE, seed = 1, num.threads = 0)
  p_rf <- predict(rf, dte)$predictions
  trk <- a_km(train); tek <- a_km(test)
  trk$r_ols <- residuals(m); trk$r_rf <- dtr$y - rf$predictions        # residuos OOB del RF
  k1 <- krigear(trk, tek, "r_ols"); k2 <- krigear(trk, tek, "r_rf")
  out <- tibble(OLS = p_ols, RK = p_ols + k1$pred, RF = p_rf, RFRK = p_rf + k2$pred)
  if (con_var) {
    nug <- function(vm) sum(vm$psill[vm$model == "Nug"])
    out$RK_var   <- k1$var + nug(k1$vm)
    out$RFRK_var <- k2$var + nug(k2$vm)
  }
  out
}
metricas <- function(o, p) c(ME = mean(p - o), MAE = mean(abs(p - o)), RMSE = sqrt(mean((p - o)^2)),
                             NSE = 1 - sum((o - p)^2) / sum((o - mean(o))^2),
                             MAE_t_ha = mean(abs(exp(p) - exp(o))) / 1000)
