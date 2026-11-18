# =============================================================================
# 05_ESDA.R — Análisis exploratorio de datos espaciales
#   (a) Tabla 1: descriptivos ponderados (UA-campaña 2015–2025)
#   (b) Moran global del rendimiento y covariables (distritos, pool 2015–2025),
#       sensibilidad a W y estabilidad por año
#   (c) LISA con permutación condicional y corrección FDR (Benjamini–Hochberg)
#   (d) Variograma empírico y ajustado a nivel de UA (2023–2024): rango -> bloques CV
# Salidas: output/tablas/T1–T3, output/figuras/F1–F3, data/interim/W_distritos.rds,
#          data/interim/variograma_ua.rds
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(spdep); library(gstat); library(ggplot2) })
set.seed(2025)
NSIM  <- 9999
N_MIN <- 15          # UA mínimas por distrito (pool plurianual)

papa <- readRDS(file.path(dir_interim, "papa_ua_cov.rds")) |> filter(ANIO >= 2015)
dsf  <- st_read(file.path(dir_interim, "distritos_papa.gpkg"), quiet = TRUE) |>
  filter(n_ua >= N_MIN) |> st_transform(CRS_M)
peru <- leer_distritos() |> st_transform(CRS_M) |> st_union()
message("Distritos analizados: ", nrow(dsf))

# ---- (a) Tabla 1: descriptivos ---------------------------------------------------
vars_t1 <- c(rend_t = "Rendimiento (t/ha)", y = "ln rendimiento (kg/ha)",
             riego = "Riego (0/1)", riego_tec = "Riego tecnificado (0/1)",
             semilla_cert = "Semilla certificada (0/1)", fert = "Fertilizante químico (0/1)",
             abono = "Abono orgánico (0/1)", plag = "Plaguicidas (0/1)",
             asist = "Asistencia técnica (0/1)", credito = "Crédito (0/1)",
             seguro = "Seguro agrario (0/1)", asoc = "Asociatividad (0/1)",
             sup_cos = "Superficie cosechada (ha)",
             ppt_d = "Precipitación de campaña (mm)", tmin_d = "Tmin media de campaña (°C)",
             tmax_d = "Tmax media de campaña (°C)", def_d = "Déficit hídrico (mm)",
             meses_helada_d = "Meses con Tmin < 0 °C", gdd_d = "Grados-día (base 4 °C)",
             alt_pt = "Altitud SRTM (m, UA con coordenadas)", pend_pt = "Pendiente (°)",
             tviaje_pt = "Tiempo a ciudad ≥ 20 mil hab. (min)")
wmean <- function(x, w) { ok <- !is.na(x) & !is.na(w); sum(w[ok] * x[ok]) / sum(w[ok]) }
wsd   <- function(x, w) { ok <- !is.na(x) & !is.na(w); m <- wmean(x, w)
                          sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }
t1 <- papa |> mutate(rend_t = rend / 1000) |>
  (\(d) tibble(variable = names(vars_t1), descripcion = unname(vars_t1),
               n = vapply(names(vars_t1), function(v) sum(!is.na(d[[v]])), 0),
               media_pond = vapply(names(vars_t1), function(v) wmean(d[[v]], d$FACTOR), 0),
               de_pond = vapply(names(vars_t1), function(v) wsd(d[[v]], d$FACTOR), 0),
               p25 = vapply(names(vars_t1), function(v) quantile(d[[v]], .25, na.rm = TRUE), 0),
               mediana = vapply(names(vars_t1), function(v) median(d[[v]], na.rm = TRUE), 0),
               p75 = vapply(names(vars_t1), function(v) quantile(d[[v]], .75, na.rm = TRUE), 0)))()
write_csv(t1, file.path(dir_tab, "T1_descriptivos.csv"))
print(t1, n = Inf, width = Inf)

# ---- (b) Matrices de pesos -----------------------------------------------------------
co <- st_coordinates(st_point_on_surface(st_geometry(dsf)))
W_knn <- function(k) nb2listw(make.sym.nb(knn2nb(knearneigh(co, k = k))), style = "W")
nb_q <- poly2nb(dsf, queen = TRUE)
islas <- which(card(nb_q) == 0)
if (length(islas)) nb_q <- make.sym.nb(union.nb(nb_q, knn2nb(knearneigh(co, k = 1))))
W <- list(knn6 = W_knn(6), knn10 = W_knn(10), reina = nb2listw(nb_q, style = "W"))
message("Contigüidad reina: ", length(islas), " islas unidas a su vecino más cercano; ",
        n.comp.nb(nb_q)$nc, " componentes")
saveRDS(list(W = W, ubigeo = dsf$ubigeo, coords = co), file.path(dir_interim, "W_distritos.rds"))

# Moran global
moran_tab <- function(x, lw, etiqueta) {
  m <- moran.mc(x, lw, nsim = NSIM, na.action = na.exclude, zero.policy = TRUE)
  tibble(variable = etiqueta, I = unname(m$statistic), E_I = -1 / (sum(!is.na(x)) - 1),
         p_sim = m$p.value)
}
vars_moran <- c("y", "riego", "semilla_cert", "fert", "asist", "credito", "asoc",
                "ppt_d", "tmin_d", "def_d", "meses_helada_d", "alt_d", "pend_d", "tviaje_d")
t2 <- bind_rows(lapply(vars_moran, function(v) moran_tab(dsf[[v]], W$knn6, v)))
t2_sens <- bind_rows(lapply(names(W), function(w) moran_tab(dsf$y, W[[w]], paste0("y [", w, "]"))))

# Estabilidad por año: medias distrito-año (>= 5 UA) con kNN-6 propio de cada año
por_anio <- papa |> filter(!is.na(ubigeo_mapa)) |>
  mutate(w = ifelse(is.na(FACTOR), 0, FACTOR)) |>
  group_by(ANIO, ubigeo = ubigeo_mapa) |>
  summarise(n = n(), y = sum(w * y) / sum(w), .groups = "drop") |>
  filter(n >= 5, is.finite(y))
cent <- st_point_on_surface(st_transform(leer_distritos(), CRS_M))
t2_anio <- bind_rows(lapply(sort(unique(por_anio$ANIO)), function(a) {
  d <- inner_join(cent, filter(por_anio, ANIO == a), by = "ubigeo")
  lw <- nb2listw(make.sym.nb(knn2nb(knearneigh(st_coordinates(d), k = 6))), style = "W")
  moran_tab(d$y, lw, paste0("y ", a)) |> mutate(n_distritos = nrow(d))
}))
t2_all <- bind_rows(mutate(t2, bloque = "Covariables (kNN-6)"),
                    mutate(t2_sens, bloque = "Sensibilidad a W"),
                    mutate(t2_anio, bloque = "Por año (kNN-6)"))
write_csv(t2_all, file.path(dir_tab, "T2_moran_global.csv"))
print(t2_all, n = Inf)

# ---- (c) LISA + FDR --------------------------------------------------------------------
lisa_clusters <- function(x, lw, alpha = 0.05) {
  li <- localmoran_perm(x, lw, nsim = NSIM)
  z  <- as.numeric(scale(x)); lz <- lag.listw(lw, z)
  p  <- li[, "Pr(z != E(Ii)) Sim"]
  q  <- p.adjust(p, "fdr")
  lab <- function(pp) factor(ifelse(pp > alpha, "No significativo",
                     ifelse(z > 0 & lz > 0, "Alto-Alto", ifelse(z < 0 & lz < 0, "Bajo-Bajo",
                     ifelse(z > 0, "Alto-Bajo", "Bajo-Alto")))),
                     levels = c("Alto-Alto", "Bajo-Bajo", "Alto-Bajo", "Bajo-Alto", "No significativo"))
  tibble(Ii = li[, "Ii"], p_sim = p, q_fdr = q, cluster_sin_fdr = lab(p), cluster = lab(q))
}
lisa <- lisa_clusters(dsf$y, W$knn6)
dsf  <- bind_cols(dsf, lisa)
t3 <- dsf |> st_drop_geometry() |> count(departamento, cluster) |>
  tidyr::pivot_wider(names_from = cluster, values_from = n, values_fill = 0)
write_csv(t3, file.path(dir_tab, "T3_LISA_por_departamento.csv"))
message("LISA (FDR 5 %): ", paste(names(table(dsf$cluster)), table(dsf$cluster), sep = "=", collapse = " | "),
        "  [sin FDR: ", sum(dsf$cluster_sin_fdr != "No significativo"), " signif.]")
st_write(dsf, file.path(dir_interim, "distritos_papa_lisa.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- Figuras 1 y 2 -------------------------------------------------------------------------
tema <- theme_void(base_size = 10) + theme(legend.position = c(0.2, 0.25),
                                          plot.title = element_text(face = "bold"))
dsf$rend_geo <- exp(dsf$y) / 1000
f1 <- ggplot() + geom_sf(data = peru, fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = dsf, aes(fill = cut(rend_geo, quantile(rend_geo, 0:5 / 5), include.lowest = TRUE)),
          colour = NA) +
  scale_fill_viridis_d("t/ha (quintiles)", option = "D") +
  labs(title = "Rendimiento de papa por distrito, 2015–2025",
       caption = sprintf("Media geométrica ponderada por FACTOR; %d distritos con ≥ %d UA. Fuente: ENA-INEI.",
                         nrow(dsf), N_MIN)) + tema
ggsave(file.path(dir_fig, "F1_rendimiento_distrital.png"), f1, width = 6, height = 7.5, dpi = 300, bg = "white")

col_lisa <- c("Alto-Alto" = "#d7191c", "Bajo-Bajo" = "#2c7bb6", "Alto-Bajo" = "#fdae61",
              "Bajo-Alto" = "#abd9e9", "No significativo" = "grey85")
f2 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = dsf, aes(fill = cluster), colour = NA) +
  scale_fill_manual("LISA (q-FDR < 0,05)", values = col_lisa, drop = FALSE) +
  labs(title = "Clústeres locales del ln-rendimiento de papa",
       caption = sprintf("Moran local, %d permutaciones, W = kNN-6; corrección de Benjamini–Hochberg.", NSIM)) + tema
ggsave(file.path(dir_fig, "F2_LISA_FDR.png"), f2, width = 6, height = 7.5, dpi = 300, bg = "white")

# ---- (d) Variograma a nivel de UA (2023–2024 = entrenamiento) -------------------------------
cov_pt <- c("alt_pt", "pend_pt", "tviaje_pt", "ppt_pt", "tmin_pt", "def_pt", "meses_helada_pt")
up <- papa |> filter(ANIO %in% 2023:2024, if_all(all_of(c("y", cov_pt)), ~ !is.na(.x))) |>
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = CRS_GEO) |> st_transform(CRS_M)
xy <- st_coordinates(up) / 1000                          # km, plano (sin CRS)
up_df <- bind_cols(st_drop_geometry(up), X_km = xy[, 1], Y_km = xy[, 2])
up_df$res <- residuals(lm(reformulate(cov_pt, "y"), data = up_df))
up_df <- st_as_sf(up_df, coords = c("X_km", "Y_km"))
vg_y   <- variogram(y ~ 1,   up_df, cutoff = 300, width = 10)
vg_res <- variogram(res ~ 1, up_df, cutoff = 300, width = 10)
ajustar <- function(vg) fit.variogram(vg, vgm(c("Exp", "Sph", "Mat")))
vm_y <- ajustar(vg_y); vm_res <- ajustar(vg_res)
rango_practico <- function(vm) {
  m <- as.character(vm$model[2]); r <- vm$range[2]
  switch(m, Exp = 3 * r, Sph = r, Mat = 4 * r, r)   # Mat (kappa~0.5–1): aproximación
}
t_vg <- tibble(variable = c("y", "residuo OLS"),
               modelo = c(as.character(vm_y$model[2]), as.character(vm_res$model[2])),
               pepita = c(vm_y$psill[1], vm_res$psill[1]),
               meseta_parcial = c(vm_y$psill[2], vm_res$psill[2]),
               rango_km = c(vm_y$range[2], vm_res$range[2]),
               rango_practico_km = c(rango_practico(vm_y), rango_practico(vm_res))) |>
  mutate(ratio_pepita = pepita / (pepita + meseta_parcial))
write_csv(t_vg, file.path(dir_tab, "T3b_variograma_UA.csv"))
print(t_vg)
saveRDS(list(vg_y = vg_y, vg_res = vg_res, vm_y = vm_y, vm_res = vm_res, cov_pt = cov_pt,
             rango_bloque_km = rango_practico(vm_res)), file.path(dir_interim, "variograma_ua.rds"))

lin <- function(vm, tag) { d <- variogramLine(vm, maxdist = 300); d$serie <- tag; d }
pts <- bind_rows(mutate(as.data.frame(vg_y)[c("dist", "gamma", "np")], serie = "ln rendimiento"),
                 mutate(as.data.frame(vg_res)[c("dist", "gamma", "np")], serie = "Residuo OLS (ambiente)"))
f3 <- ggplot(pts, aes(dist, gamma, colour = serie)) + geom_point(aes(size = np), alpha = .7) +
  geom_line(data = bind_rows(lin(vm_y, "ln rendimiento"), lin(vm_res, "Residuo OLS (ambiente)"))) +
  scale_size_continuous("Pares", range = c(0.8, 3)) + scale_colour_manual(NULL, values = c("#1b9e77", "#d95f02")) +
  labs(x = "Distancia (km)", y = "Semivarianza", title = "Variograma empírico y ajustado (UA 2023–2024)",
       caption = sprintf("n = %d UA; coordenadas de segmento. Rango práctico del residuo ≈ %.0f km.",
                         nrow(up), rango_practico(vm_res))) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(dir_fig, "F3_variograma_UA.png"), f3, width = 7, height = 4.5, dpi = 300, bg = "white")
message("ESDA listo.")
