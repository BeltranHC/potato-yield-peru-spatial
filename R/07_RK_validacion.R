# =============================================================================
# 07_RK_validacion.R — Etapa 3: predicción a nivel de UA y validación honesta
# Entrenamiento: UA 2023–2024 con coordenadas | Prueba externa: UA 2025 (año nunca visto)
# Modelos:
#   OLS, RK (OLS + kriging de residuos), RF (ranger) y RFRK (RF + kriging de residuos OOB)
#   con dos conjuntos de predictores:
#     A = solo ambientales (exhaustivos -> mapeables): altitud, pendiente, ln tiempo a ciudad,
#         ppt, tmin, déficit, meses con helada (ventana de campaña, buffer 1 km)
#     B = A + manejo de la ENA (no mapeables; predicción para UA encuestadas)
# Kriging: la pepita (68 % de la varianza) se modela como error de medición ("Err"):
#   las UA de un mismo segmento comparten coordenada, y una UA nueva en el mismo
#   punto no debe recibir el valor exacto de sus vecinas.
# Validación cruzada (k = 10) dentro de 2023–2024:
#   aleatoria | bloques espaciales (tamaño = rango del variograma, blockCV) |
#   kNNDM (CAST, ajusta la distribución de distancias al dominio de predicción) | por año
# Métricas en escala log: ME, MAE, RMSE, NSE; y MAE en t/ha.
# Salidas: output/tablas/T10, T10b; output/figuras/F6 (distancias), F7 (mapa)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(terra); library(gstat); library(ranger)
  library(CAST); library(blockCV); library(exactextractr); library(ggplot2); library(patchwork) })
set.seed(2025)
source("R/funciones_prediccion.R")      # K, NTREE, NMAX, A, B, preparar_puntos(), krigear(), ...

papa <- readRDS(file.path(dir_interim, "papa_ua_cov.rds"))
pt <- preparar_puntos(papa)
tr <- pt[pt$ANIO %in% 2023:2024, ]; te <- pt[pt$ANIO == 2025, ]
message("Entrenamiento 2023–2024: ", nrow(tr), " UA | Prueba 2025: ", nrow(te), " UA")

# ---- Dominio de predicción: distritos con papa ---------------------------------------------------
dominio <- st_read(file.path(dir_interim, "distritos_papa.gpkg"), quiet = TRUE) |>
  st_transform(CRS_M) |> st_union() |> st_make_valid() |> st_simplify(dTolerance = 500)

# ---- Pliegues ---------------------------------------------------------------------------------------
f_folds <- file.path(dir_interim, "folds_ua_2023_2024.rds")
if (!file.exists(f_folds)) {
  folds <- list(aleatoria = sample(rep(1:K, length.out = nrow(tr))))
  rango_m <- vgi$rango_bloque_km * 1000
  folds$bloques <- cv_spatial(x = tr, k = K, size = rango_m, selection = "random", iteration = 50,
                              seed = 2025, progress = FALSE, report = FALSE, plot = FALSE)$folds_ids
  # kNNDM sobre ubicaciones únicas (UA del mismo segmento comparten coordenada y pliegue)
  loc <- paste(round(tr$X, 3), round(tr$Y, 3))
  u   <- tr[!duplicated(loc), ]
  kn  <- knndm(tpoints = u, modeldomain = dominio, k = K, clustering = "hierarchical", samplesize = 2000)
  folds$kNNDM <- kn$clusters[match(loc, loc[!duplicated(loc)])]
  message("kNNDM: ", nrow(u), " ubicaciones únicas; W = ", round(kn$W, 4))
  folds$anio <- as.integer(factor(tr$ANIO))
  saveRDS(folds, f_folds)
}
folds <- readRDS(f_folds)
print(sapply(folds, function(f) length(unique(f))))

# ---- Validación cruzada ------------------------------------------------------------------------------
res_cv <- list()
for (esq in names(folds)) for (set in c("A", "B")) {
  t0 <- Sys.time(); fid <- folds[[esq]]; vars <- get(set)
  pr <- matrix(NA_real_, nrow(tr), 4, dimnames = list(NULL, c("OLS", "RK", "RF", "RFRK")))
  for (f in sort(unique(fid))) {
    k <- fid == f
    pr[k, ] <- as.matrix(ajustar_predecir(tr[!k, ], tr[k, ], vars))
  }
  res_cv[[paste(esq, set)]] <- bind_rows(lapply(colnames(pr), function(mn)
    as_tibble(as.list(metricas(tr$y, pr[, mn]))) |> mutate(validacion = esq, predictores = set, modelo = mn)))
  message(sprintf("CV %-9s | %s | %.1f min", esq, set, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
# Prueba externa: entrenar con 2023–2024, predecir 2025
for (set in c("A", "B")) {
  p25 <- ajustar_predecir(tr, te, get(set))
  res_cv[[paste("2025", set)]] <- bind_rows(lapply(names(p25), function(mn)
    as_tibble(as.list(metricas(te$y, p25[[mn]]))) |>
      mutate(validacion = "prueba 2025", predictores = set, modelo = mn)))
}
t10 <- bind_rows(res_cv) |> relocate(validacion, predictores, modelo) |>
  mutate(validacion = factor(validacion, c("aleatoria", "bloques", "kNNDM", "anio", "prueba 2025"))) |>
  arrange(predictores, modelo, validacion)
write_csv(t10, file.path(dir_tab, "T10_validacion_RK.csv"))
print(t10, n = Inf, width = Inf)

# Inflación: cuánto sobreestima la CV aleatoria el NSE respecto a la prueba 2025
t10b <- t10 |> select(validacion, predictores, modelo, NSE, RMSE) |>
  tidyr::pivot_wider(names_from = validacion, values_from = c(NSE, RMSE)) |>
  mutate(inflacion_NSE_aleatoria = `NSE_aleatoria` - `NSE_prueba 2025`,
         inflacion_NSE_kNNDM = `NSE_kNNDM` - `NSE_prueba 2025`)
write_csv(t10b, file.path(dir_tab, "T10b_inflacion_validacion.csv"))
print(t10b, width = Inf)

# ---- Figura 6: distancias de vecino más cercano (CAST::geodist) ----------------------------------
gd <- bind_rows(lapply(c("aleatoria", "bloques", "kNNDM"), function(esq) {
  g <- tryCatch(geodist(tr, modeldomain = dominio, CVtest = folds[[esq]], testdata = te,
                        samplesize = 2000), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  as.data.frame(g) |> mutate(esquema = esq)
}))
if (nrow(gd)) {
  gd <- gd |> mutate(what = as.character(what),
    serie = case_when(what == "sample-to-sample" ~ "muestra–muestra",
                      what == "prediction-to-sample" ~ "dominio–muestra",
                      what == "test-to-sample" ~ "2025–entrenamiento",
                      grepl("CV", what) ~ paste0("CV ", esquema), TRUE ~ what)) |>
    distinct(serie, dist, .keep_all = TRUE)
  f6 <- ggplot(gd, aes(dist / 1000, colour = serie)) + geom_density(linewidth = .8) +
    scale_x_continuous("Distancia al vecino más cercano en entrenamiento (km)", trans = "log1p",
                       breaks = c(0, 1, 5, 10, 25, 50, 100, 250)) +
    labs(y = "Densidad", colour = NULL,
         title = "¿Imita la validación cruzada la distancia real de predicción?") +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(dir_fig, "F6_distancias_CV.png"), f6, width = 8, height = 4.8, dpi = 300, bg = "white")
}

# ---- Figura 7: mapa de predicción RK-A (campaña oct 2024 – abr 2025) -------------------------------
grid <- st_make_grid(dominio, cellsize = 5000, what = "centers") |> st_sf() |> st_filter(st_sf(dominio))
buf  <- st_buffer(grid, 1000) |> st_transform(CRS_GEO)
dem  <- rast(file.path(dir_ext, "dem", "srtm90_andes.tif")); pend <- rast(file.path(dir_ext, "dem", "pendiente90_andes.tif"))
acc  <- rast(list.files(file.path(dir_ext, "acceso"), "travel_time.*[.]tif$", recursive = TRUE, full.names = TRUE)[1])
g <- tibble(alt_pt = exact_extract(dem, buf, "mean", progress = FALSE),
            pend_pt = exact_extract(pend, buf, "mean", progress = FALSE),
            ln_tviaje = log1p(exact_extract(acc, buf, "mean", progress = FALSE)))
camp <- function(v, fun) {
  r <- rast(file.path(dir_ext, "terraclimate", sprintf("tc_%s_%d.nc", v, 2024:2025)))
  r <- r[[10:16]]                                   # oct 2024 ... abr 2025
  m <- as.matrix(exact_extract(r, buf, "mean", progress = FALSE)); apply(m, 1, fun)
}
g$ppt_pt <- camp("ppt", sum); g$tmin_pt <- camp("tmin", mean); g$def_pt <- camp("def", sum)
g$meses_helada_pt <- camp("tmin", function(x) sum(x < 0))
xyg <- st_coordinates(grid) / 1000; g$X <- xyg[, 1]; g$Y <- xyg[, 2]
ok <- complete.cases(g)
todo <- rbind(tr, te)
m_all <- lm(reformulate(A, "y"), st_drop_geometry(todo))
tk <- a_km(todo); tk$r <- residuals(m_all)
gk <- st_as_sf(g[ok, ], coords = c("X", "Y"), remove = FALSE)
kr <- krigear(tk, gk, "r")
mapa <- st_sf(geometry = st_geometry(grid)[ok], rend_t = exp(predict(m_all, g[ok, ]) + kr$pred) / 1000,
              sd_krig = sqrt(kr$var))
# Área de aplicabilidad (Meyer & Pebesma 2021): DI ponderado por importancia del RF,
# umbral calibrado con los pliegues kNNDM de 2023–2024. Fuera del AOA = extrapolación.
rf_imp <- ranger(reformulate(A, "y"), data = st_drop_geometry(tr)[c("y", A)], num.trees = NTREE,
                 importance = "permutation", seed = 1)
aoa_r <- aoa(newdata = as.data.frame(g[ok, A]), train = as.data.frame(st_drop_geometry(tr)[A]),
             variables = A, weight = as.data.frame(t(pmax(rf_imp$variable.importance, 0))),
             CVtest = folds$kNNDM, verbose = FALSE)
mapa$DI <- aoa_r$DI; mapa$AOA <- aoa_r$AOA == 1
message("Celdas dentro del AOA: ", round(100 * mean(mapa$AOA), 1), " %")
saveRDS(mapa, file.path(dir_interim, "mapa_RK_2025.rds"))
mapa_in <- mapa[mapa$AOA, ]; mapa_out <- mapa[!mapa$AOA, ]
peru <- leer_distritos() |> st_transform(CRS_M) |> st_union()
base <- function() list(geom_sf(data = peru, fill = "white", colour = "grey60", linewidth = .2),
                        theme_void(base_size = 9))
fuera <- geom_sf(data = mapa_out, colour = "grey80", size = .25, shape = 15)
f7 <- (ggplot() + base() + fuera + geom_sf(data = mapa_in, aes(colour = rend_t), size = .25, shape = 15) +
         scale_colour_viridis_c("t/ha", option = "D") + labs(title = "A. Rendimiento predicho (RK, ambiente)")) +
      (ggplot() + base() + fuera + geom_sf(data = mapa_in, aes(colour = sd_krig), size = .25, shape = 15) +
         scale_colour_viridis_c("DE kriging\n(log)", option = "magma") + labs(title = "B. Incertidumbre del kriging")) +
      plot_annotation(caption = paste0("Rejilla de 5 km sobre distritos con papa; campaña oct 2024 – abr 2025. ",
        "Gris: fuera del área de aplicabilidad (AOA; Meyer & Pebesma, 2021)."))
ggsave(file.path(dir_fig, "F7_mapa_RK_2025_AOA.png"), f7, width = 10, height = 6.5, dpi = 300, bg = "white")
message("Etapa 3 lista.")
