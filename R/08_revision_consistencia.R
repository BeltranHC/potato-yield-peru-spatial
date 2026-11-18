# =============================================================================
# 08_revision_consistencia.R — respuestas a revisores (comentarios 1 y 5)
#  (1) Consistencia entre escalas: el efecto de cada práctica de manejo estimado
#      (a) entre distritos (OLS, SDM directo, mediana MGWR) en unidades de "0 → 100 %
#      de adopción", y (b) dentro del distrito a nivel de UA (EF distrito + año + mes de
#      siembra, clima de la campaña; EE agrupados por distrito). Si los signos coinciden,
#      la lectura distrital es confiable; si no, hay falacia ecológica.
#  (5) Robustez del SDM distrital: umbral de UA, marco muestral antiguo vs. nuevo,
#      sin 2021, umbral de 50 UA y panel distrito-año con EF de año.
# Salidas: output/tablas/T11_consistencia_escalas.csv, T12_robustez_SDM.csv,
#          output/figuras/F8_consistencia_escalas.png
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(spdep); library(spatialreg); library(fixest)
                                 library(ggplot2) })
set.seed(2025)

papa <- readRDS(file.path(dir_interim, "papa_ua_cov.rds")) |> filter(ANIO >= 2015, !is.na(ubigeo_mapa))
md   <- readRDS(file.path(dir_interim, "modelos_distritales.rds"))
dsf  <- st_read(file.path(dir_interim, "distritos_papa_lisa.gpkg"), quiet = TRUE) |>
  mutate(ln_tviaje = log1p(tviaje_d))
X <- md$X; etq <- md$etiquetas
manejo <- c("riego", "semilla_cert", "fert", "asist", "credito", "asoc", "ln_sup")
sd_x <- sapply(X, function(v) sd(dsf[[v]], na.rm = TRUE))      # DE distrital (para des-estandarizar)
pct <- function(b) 100 * (exp(b) - 1)

# ---- (1) Consistencia entre escalas --------------------------------------------------------------
# (a) Distrital, en unidades naturales: β_z / DE  =  efecto de 0 → 1 en la proporción de adopción
ols_ct <- summary(md$ols)$coefficients
imp <- read_csv(file.path(dir_tab, "T7_impactos_SDM_SDEM.csv"), show_col_types = FALSE) |>
  filter(modelo == "SDM", efecto == "direct")
f_mg <- file.path(dir_tab, "T16_MGWR_inferencia.csv")        # mgwr (inferencia exacta); si no, GWmodel
if (!file.exists(f_mg)) f_mg <- file.path(dir_tab, "T8_MGWR_resumen.csv")
mg <- read_csv(f_mg, show_col_types = FALSE)
dist_tab <- tibble(variable = manejo,
  ols_b = ols_ct[manejo, 1] / sd_x[manejo], ols_p = ols_ct[manejo, 4],
  sdm_b = imp$estimacion[match(manejo, imp$variable)] / sd_x[manejo], sdm_p = imp$p[match(manejo, imp$variable)],
  mgwr_med = mg$mediana[match(manejo, mg$variable)] / sd_x[manejo])

# (b) UA dentro del distrito (todas las prácticas a la vez)
u <- papa |>
  mutate(mes_siembra = ifelse(ventana_declarada, as.integer(format(ini, "%m")), 0L),
         cred_na = as.integer(is.na(credito)),
         across(c(credito, semilla_cert, fert, abono, plag, asist, asoc), ~ ifelse(is.na(.x), 0, .x)),
         z_def = as.numeric(scale(def_d)), z_tmin = as.numeric(scale(tmin_d)), z_hel = as.numeric(scale(meses_helada_d)),
         segmento = sub("^([0-9]+_[^_]+)_.*$", "\\1", id_ua))
rhs <- "riego + semilla_cert + fert + asist + credito + asoc + ln_sup + abono + plag + cred_na"
m_ua  <- feols(as.formula(paste("y ~", rhs, "+ z_def + z_tmin + z_hel | ubigeo_mapa + ANIO + mes_siembra")),
               u, cluster = ~ubigeo_mapa)
m_seg <- feols(as.formula(paste("y ~", rhs, "| segmento + mes_siembra")),
               filter(u, ANIO >= 2023), cluster = ~ubigeo_mapa)
ct_ua <- coeftable(m_ua); ct_seg <- coeftable(m_seg)
t11 <- dist_tab |>
  mutate(ua_b = ct_ua[manejo, 1], ua_p = ct_ua[manejo, 4],
         seg_b = ct_seg[manejo, 1], seg_p = ct_seg[manejo, 4],
         etiqueta = etq[variable],
         dist_p = pmin(ols_p, sdm_p),               # significativo en OLS o en SDM
         signo_coincide = sign(sdm_b) == sign(ua_b),
         lectura = case_when(
           dist_p < 0.05 & ua_p < 0.05 & sign(sdm_b) == sign(ua_b) ~ "Consistente",
           dist_p < 0.05 & ua_p < 0.05 & sign(sdm_b) != sign(ua_b) ~ "Falacia ecológica (signo opuesto)",
           dist_p < 0.05 & ua_p >= 0.05 ~ "Solo entre distritos (composición/contexto)",
           dist_p >= 0.05 & ua_p < 0.05 ~ "Solo dentro del distrito",
           TRUE ~ "No significativo en ninguna escala")) |>
  relocate(etiqueta)
write_csv(t11, file.path(dir_tab, "T11_consistencia_escalas.csv"))
print(t11 |> mutate(across(c(ols_b, sdm_b, mgwr_med, ua_b, seg_b), ~ round(pct(.x), 1))), width = Inf)

# Figura 8: distrital (SDM directo) vs. UA (EF distrito), en %
fd <- bind_rows(
  tibble(etiqueta = etq[manejo], escala = "Entre distritos (SDM, efecto directo)",
         b = t11$sdm_b, ee = NA_real_),
  tibble(etiqueta = etq[manejo], escala = "Dentro del distrito (UA, EF distrito)",
         b = ct_ua[manejo, 1], ee = ct_ua[manejo, 2]),
  tibble(etiqueta = etq[manejo], escala = "Dentro del segmento (UA 2023–2025, ~1 km)",
         b = ct_seg[manejo, 1], ee = ct_seg[manejo, 2]))
se_sdm <- imp$estimacion[match(manejo, imp$variable)] / imp$z[match(manejo, imp$variable)] / sd_x[manejo]
fd$ee[fd$escala == "Entre distritos (SDM, efecto directo)"] <- se_sdm
fd$etiqueta <- factor(fd$etiqueta, levels = rev(etq[manejo]))
f8 <- ggplot(fd, aes(pct(b), etiqueta, colour = escala)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = pct(b - 1.96 * ee), xmax = pct(b + 1.96 * ee)), width = .25,
                position = position_dodge(width = .6), orientation = "y") +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  scale_colour_manual(NULL, values = c("#b2182b", "#2166ac", "#1b7837")) +
  labs(x = "Efecto sobre el rendimiento (%; adopción 0 → 1; IC 95 %)", y = NULL,
       title = "Consistencia del efecto del manejo entre escalas") +
  coord_cartesian(xlim = c(-60, 120)) + theme_minimal(base_size = 9) + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
ggsave(file.path(dir_fig, "F8_consistencia_escalas.png"), f8, width = 7.5, height = 5, dpi = 300, bg = "white")

# ---- (5) Robustez del SDM distrital -----------------------------------------------------------------
estat <- st_drop_geometry(dsf)[c("ubigeo", "alt_d", "pend_d", "ln_tviaje")]
cent  <- st_point_on_surface(st_transform(dsf["ubigeo"], CRS_M))
covs_ua <- c("y", "riego", "semilla_cert", "fert", "asist", "credito", "asoc", "ln_sup",
             "def_d", "tmin_d", "meses_helada_d")
agregar <- function(d, por = "ubigeo") {
  d |> mutate(w = ifelse(is.na(FACTOR), 0, FACTOR), ubigeo = ubigeo_mapa) |>
    group_by(across(all_of(por))) |>
    summarise(n_ua = n(), across(all_of(covs_ua), ~ { ok <- !is.na(.x) & w > 0
      if (any(ok)) weighted.mean(.x[ok], w[ok]) else NA_real_ }), .groups = "drop") |>
    left_join(estat, by = "ubigeo")
}
X_nat <- X                                   # unidades naturales (proporciones 0–1, mm, °C...)
sdm_variante <- function(dd, nombre, n_min = 15, W = NULL, extra = NULL) {
  dd <- dd |> filter(n_ua >= n_min) |> filter(if_all(all_of(c("y", X_nat)), ~ is.finite(.x)))
  if (is.null(W)) {
    co <- st_coordinates(cent)[match(dd$ubigeo, cent$ubigeo), ]
    W  <- nb2listw(make.sym.nb(knn2nb(knearneigh(co, k = 6))), style = "W")
  }
  f <- reformulate(c(X_nat, extra), "y")
  m <- lagsarlm(f, dd, W, Durbin = reformulate(X_nat))
  tr <- trW(as(W, "CsparseMatrix"), type = "MC")
  s <- summary(impacts(m, tr = tr, R = 1000), zstats = TRUE, short = TRUE)
  bind_rows(lapply(1:2, function(j) {
    k <- c("direct", "indirect")[j]; ix <- match(X_nat, names(s$res[[k]]))
    e <- s$res[[k]][ix]; se <- s$semat[ix, j]
    tibble(variante = nombre, n = nrow(dd), rho = m$rho, efecto = k, variable = X_nat,
           b = unname(e), p = 2 * pnorm(-abs(unname(e / se))))
  }))
}
base <- agregar(papa)
variantes <- list(
  sdm_variante(base, "Base (≥ 15 UA, 2015–2025)"),
  sdm_variante(base, "≥ 30 UA", n_min = 30),
  sdm_variante(agregar(filter(papa, ANIO <= 2022)), "Marco antiguo (2015–2022)"),
  sdm_variante(agregar(filter(papa, ANIO >= 2023)), "Marco nuevo (2023–2025)"),
  sdm_variante(agregar(filter(papa, ANIO != 2021)), "Sin 2021"),
  sdm_variante(base, "≥ 50 UA", n_min = 50))
# Panel distrito-año (≥ 5 UA), W bloque-diagonal (kNN-6 dentro de cada año) y EF de año
pan <- agregar(papa, por = c("ubigeo", "ANIO")) |> filter(n_ua >= 5) |>
  filter(if_all(all_of(c("y", X_nat)), ~ is.finite(.x))) |> arrange(ANIO, ubigeo)
nb_list <- lapply(split(seq_len(nrow(pan)), pan$ANIO), function(ix) {
  co <- st_coordinates(cent)[match(pan$ubigeo[ix], cent$ubigeo), ]
  make.sym.nb(knn2nb(knearneigh(co, k = 6)))
})
off <- cumsum(c(0, head(sapply(nb_list, length), -1)))
nb_pan <- structure(unlist(lapply(seq_along(nb_list), function(i) lapply(nb_list[[i]], function(v) as.integer(v + off[i]))),
                           recursive = FALSE), class = "nb", region.id = as.character(seq_len(nrow(pan))))
pan$f_anio <- factor(pan$ANIO)
variantes[[length(variantes) + 1]] <- sdm_variante(pan, "Panel distrito-año + EF año", n_min = 5,
                                                   W = nb2listw(nb_pan, style = "W"), extra = "f_anio")
t12 <- bind_rows(variantes) |> mutate(etiqueta = etq[variable])
write_csv(t12, file.path(dir_tab, "T12_robustez_SDM.csv"))
print(t12 |> filter(variable %in% c("riego", "semilla_cert", "fert", "credito", "ln_sup", "tmin_d")) |>
        mutate(v = sprintf("%.1f%s", 100 * b, ifelse(p < .05, "*", ""))) |>
        select(variante, n, efecto, variable, v) |>
        tidyr::pivot_wider(names_from = c(variable), values_from = v), n = Inf, width = Inf)
message("08 listo.")
