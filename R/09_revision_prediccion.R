# =============================================================================
# 09_revision_prediccion.R — respuestas a revisores (comentarios 2 y 4)
#  (2) ¿De dónde viene la inflación de la CV aleatoria? Se agregan dos CV agrupadas:
#        "segmento-año": deja fuera segmentos completos de una campaña (sin vecinas
#                        del mismo segmento y año en entrenamiento)
#        "ubicación":    deja fuera todas las campañas de una misma coordenada
#      Si la inflación desaparece con la agrupación, la fuga se debe a UA que
#      comparten segmento (y campaña), no a la proximidad geográfica en general.
#  (4) Calidad de la incertidumbre y de la retrotransformación:
#        - cobertura del intervalo de predicción al 90 % y MSDR del kriging (RK, RFRK)
#        - corrección "smearing" de Duan para las métricas en t/ha
#        - bootstrap por conglomerado (segmento) de las diferencias de RMSE en 2025
# Salidas: output/tablas/T13_cv_agrupada.csv, T14_calibracion.csv, T15_bootstrap_2025.csv
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(gstat); library(ranger) })
set.seed(2025)
source("R/funciones_prediccion.R")

papa <- readRDS(file.path(dir_interim, "papa_ua_cov.rds"))
pt <- preparar_puntos(papa)
tr <- pt[pt$ANIO %in% 2023:2024, ]; te <- pt[pt$ANIO == 2025, ]
folds <- readRDS(file.path(dir_interim, "folds_ua_2023_2024.rds"))
stopifnot(length(folds$kNNDM) == nrow(tr))

agrupar <- function(g) { ug <- unique(g); f <- sample(rep(1:K, length.out = length(ug))); f[match(g, ug)] }
folds$`segmento-año` <- agrupar(tr$segmento)
folds$`ubicación`    <- agrupar(tr$ubicacion)
message("Segmentos-año: ", n_distinct(tr$segmento), " | ubicaciones: ", n_distinct(tr$ubicacion))

cv_pred <- function(fid, vars, con_var = FALSE) {
  out <- NULL
  for (f in sort(unique(fid))) {
    k <- fid == f
    p <- ajustar_predecir(tr[!k, ], tr[k, ], vars, con_var = con_var)
    if (is.null(out)) out <- as.data.frame(matrix(NA_real_, nrow(tr), ncol(p), dimnames = list(NULL, names(p))))
    out[k, ] <- as.data.frame(p)
  }
  as_tibble(out)
}

# ---- (2) CV agrupada ---------------------------------------------------------------------------
res <- list(); pred_knndm <- list()
for (esq in c("aleatoria", "segmento-año", "ubicación", "kNNDM")) for (set in c("A", "B")) {
  t0 <- Sys.time()
  p <- cv_pred(folds[[esq]], get(set), con_var = (esq == "kNNDM"))
  if (esq == "kNNDM") pred_knndm[[set]] <- p
  res[[paste(esq, set)]] <- bind_rows(lapply(c("OLS", "RK", "RF", "RFRK"), function(mn)
    as_tibble(as.list(metricas(tr$y, p[[mn]]))) |> mutate(validacion = esq, predictores = set, modelo = mn)))
  message(sprintf("CV %-13s | %s | %.1f min", esq, set, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
p25 <- lapply(c(A = "A", B = "B"), function(set) ajustar_predecir(tr, te, get(set), con_var = TRUE))
for (set in c("A", "B")) res[[paste("2025", set)]] <- bind_rows(lapply(c("OLS", "RK", "RF", "RFRK"), function(mn)
  as_tibble(as.list(metricas(te$y, p25[[set]][[mn]]))) |> mutate(validacion = "prueba 2025", predictores = set, modelo = mn)))
t13 <- bind_rows(res) |> relocate(validacion, predictores, modelo)
write_csv(t13, file.path(dir_tab, "T13_cv_agrupada.csv"))
print(t13 |> select(validacion, predictores, modelo, NSE) |>
        tidyr::pivot_wider(names_from = validacion, values_from = NSE), width = Inf)

# ---- (4a) Calibración de la incertidumbre del kriging ----------------------------------------------
calib <- function(o, p, v, nombre) {
  z <- (o - p) / sqrt(v)
  tibble(conjunto = nombre, cobertura_90 = mean(abs(z) <= qnorm(.95)), cobertura_50 = mean(abs(z) <= qnorm(.75)),
         MSDR = mean(z^2), ancho_medio_90_log = mean(2 * qnorm(.95) * sqrt(v)))
}
# ---- (4b) Smearing de Duan: factor estimado con los residuos de la CV kNNDM (entrenamiento) ------
smear <- function(o, p) mean(exp(o - p))
t14 <- bind_rows(lapply(c("A", "B"), function(set) bind_rows(lapply(c("RK", "RFRK"), function(mn) {
  pk <- pred_knndm[[set]]; pp <- p25[[set]]; vname <- paste0(mn, "_var")
  s <- smear(tr$y, pk[[mn]])
  bind_rows(calib(tr$y, pk[[mn]], pk[[vname]], "CV kNNDM 2023–2024"),
            calib(te$y, pp[[mn]], pp[[vname]], "Prueba 2025")) |>
    mutate(predictores = set, modelo = mn, factor_smearing = s,
           ME_t_ha_sin = c(mean(exp(pk[[mn]]) - exp(tr$y)), mean(exp(pp[[mn]]) - exp(te$y))) / 1000,
           ME_t_ha_con = c(mean(s * exp(pk[[mn]]) - exp(tr$y)), mean(s * exp(pp[[mn]]) - exp(te$y))) / 1000,
           MAE_t_ha_sin = c(mean(abs(exp(pk[[mn]]) - exp(tr$y))), mean(abs(exp(pp[[mn]]) - exp(te$y)))) / 1000,
           MAE_t_ha_con = c(mean(abs(s * exp(pk[[mn]]) - exp(tr$y))), mean(abs(s * exp(pp[[mn]]) - exp(te$y)))) / 1000)
}))))  |> relocate(predictores, modelo, conjunto)
write_csv(t14, file.path(dir_tab, "T14_calibracion.csv"))
print(t14, width = Inf)

# ---- (4c) Bootstrap por segmento de las diferencias de RMSE en 2025 ---------------------------------
comparaciones <- list(
  `RK − OLS (A)`        = list(c("A", "RK"),   c("A", "OLS")),
  `RFRK − RF (A)`       = list(c("A", "RFRK"), c("A", "RF")),
  `RF − OLS (A)`        = list(c("A", "RF"),   c("A", "OLS")),
  `RFRK − RK (A)`       = list(c("A", "RFRK"), c("A", "RK")),
  `RFRK (B) − RFRK (A)` = list(c("B", "RFRK"), c("A", "RFRK")),
  `RFRK (B) − RK (B)`   = list(c("B", "RFRK"), c("B", "RK")))
seg <- te$segmento; useg <- unique(seg); idx <- split(seq_along(seg), seg)
rmse <- function(e) sqrt(mean(e^2))
NB <- 2000
bs <- replicate(NB, { s <- sample(useg, replace = TRUE); unlist(idx[s], use.names = FALSE) })
t15 <- bind_rows(lapply(names(comparaciones), function(nm) {
  a <- comparaciones[[nm]][[1]]; b <- comparaciones[[nm]][[2]]
  ea <- te$y - p25[[a[1]]][[a[2]]]; eb <- te$y - p25[[b[1]]][[b[2]]]
  d  <- vapply(bs, function(ix) rmse(ea[ix]) - rmse(eb[ix]), 0)
  tibble(comparacion = nm, dif_RMSE = rmse(ea) - rmse(eb), ic95_inf = quantile(d, .025), ic95_sup = quantile(d, .975),
         p_bootstrap = 2 * min(mean(d >= 0), mean(d <= 0)))
}))
write_csv(t15, file.path(dir_tab, "T15_bootstrap_2025.csv"))
print(t15, width = Inf)
message("09 listo.")
