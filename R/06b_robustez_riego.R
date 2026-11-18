# =============================================================================
# 06b_robustez_riego.R — ¿Por qué el riego sale negativo a nivel distrital?
# Hipótesis: es una falacia ecológica / selección. El riego se concentra donde
# falta agua (r riego–déficit = 0,55 entre distritos), así que marca zonas difíciles
# y no causa el bajo rendimiento.
# Estrategia:
#   (A) Distritos: interacción riego × déficit hídrico (OLS y SEM).
#   (B) UA (2015–2025): comparación DENTRO del distrito (efectos fijos), agregando
#       año, mes de siembra (campaña grande vs. chica), clima de la campaña,
#       interacciones con déficit y helada, segmento×año (2023–2025, vecinos a ~1 km),
#       pesos FACTOR, región natural y tipo de riego (tecnificado vs. gravedad).
# EE agrupados por distrito en todas las especificaciones UA.
# Salidas: output/tablas/T9_robustez_riego.csv, T9b_riego_distrital_interaccion.csv,
#          output/figuras/F5_robustez_riego.png
# =============================================================================
source("R/00_config.R")
instalar_faltantes("fixest")
suppressPackageStartupMessages({ library(fixest); library(sf); library(spdep); library(spatialreg)
                                 library(ggplot2); library(patchwork) })

# ---- (A) Distritos: interacción riego × déficit ------------------------------------------
dsf <- st_read(file.path(dir_interim, "distritos_papa_lisa.gpkg"), quiet = TRUE) |>
  mutate(ln_tviaje = log1p(tviaje_d))
md  <- readRDS(file.path(dir_interim, "modelos_distritales.rds"))
lw  <- readRDS(file.path(dir_interim, "W_distritos.rds"))$W$knn6
d   <- st_drop_geometry(dsf); d[md$X] <- scale(d[md$X])
d$riego_x_def    <- d$riego * d$def_d
d$riego_x_helada <- d$riego * d$meses_helada_d
f_int <- reformulate(c(md$X, "riego_x_def", "riego_x_helada"), "y")
ols_i <- lm(f_int, d); sem_i <- errorsarlm(f_int, d, lw)
tb <- function(m, nm) { s <- if (inherits(m, "lm")) summary(m)$coefficients else summary(m)$Coef
  tibble(modelo = nm, termino = rownames(s), b = s[, 1], ee = s[, 2], p = s[, 4]) |>
    filter(grepl("riego|def_d|helada", termino)) }
t9b <- bind_rows(tb(ols_i, "OLS distrital + interacción"), tb(sem_i, "SEM distrital + interacción"))
write_csv(t9b, file.path(dir_tab, "T9b_riego_distrital_interaccion.csv"))
print(t9b)

# ---- (B) Nivel UA ---------------------------------------------------------------------------
p <- readRDS(file.path(dir_interim, "papa_ua_cov.rds")) |>
  filter(ANIO >= 2015, !is.na(ubigeo_mapa)) |>
  mutate(distrito = ubigeo_mapa,
         segmento = sub("^([0-9]+_[^_]+)_.*$", "\\1", id_ua),      # ANIO_NSEGM (2023+) o ANIO_CONGLOMERADO
         mes_siembra = ifelse(ventana_declarada, as.integer(format(ini, "%m")), NA),
         mes_siembra = ifelse(is.na(mes_siembra), 0L, mes_siembra),  # 0 = no declarado
         cred_na = as.integer(is.na(credito)), credito = ifelse(is.na(credito), 0, credito),
         across(c(semilla_cert, fert, abono, plag, asist, asoc, seguro), ~ ifelse(is.na(.x), 0, .x)),
         riego_grav = as.integer(riego == 1 & riego_tec == 0),        # gravedad / otros
         sierra = as.integer(REGION == 2),
         z_def = as.numeric(scale(def_d)), z_hel = as.numeric(scale(meses_helada_d)),
         z_tmin = as.numeric(scale(tmin_d)))
message("UA en el análisis: ", nrow(p), " | distritos: ", n_distinct(p$distrito))

ctrl <- "semilla_cert + fert + abono + plag + asist + credito + cred_na + asoc + ln_sup"
clim <- "z_def + z_tmin + z_hel"
esp <- list(
  `1 Pooled, sin controles`                   = feols(as.formula("y ~ riego"), p, cluster = ~distrito),
  `2 + manejo + EF año`                       = feols(as.formula(paste("y ~ riego +", ctrl, "| ANIO")), p, cluster = ~distrito),
  `3 + EF distrito (dentro del distrito)`     = feols(as.formula(paste("y ~ riego +", ctrl, "| distrito + ANIO")), p, cluster = ~distrito),
  `4 + clima de campaña + EF mes de siembra`  = feols(as.formula(paste("y ~ riego +", ctrl, "+", clim, "| distrito + ANIO + mes_siembra")), p, cluster = ~distrito),
  `5 + riego × déficit, riego × helada`       = feols(as.formula(paste("y ~ riego + riego:z_def + riego:z_hel +", ctrl, "+", clim, "| distrito + ANIO + mes_siembra")), p, cluster = ~distrito),
  `6 Tecnificado vs. gravedad`                = feols(as.formula(paste("y ~ riego_tec + riego_grav +", ctrl, "+", clim, "| distrito + ANIO + mes_siembra")), p, cluster = ~distrito),
  `7 Riego × sierra`                          = feols(as.formula(paste("y ~ riego + riego:sierra +", ctrl, "+", clim, "| distrito + ANIO + mes_siembra")), p, cluster = ~distrito),
  `8 Ponderado por FACTOR`                    = feols(as.formula(paste("y ~ riego +", ctrl, "+", clim, "| distrito + ANIO + mes_siembra")), p, weights = ~FACTOR, cluster = ~distrito),
  `9 EF segmento×año (2023–2025, ~1 km)`      = feols(as.formula(paste("y ~ riego +", ctrl, "| segmento + mes_siembra")), filter(p, ANIO >= 2023), cluster = ~distrito)
)
t9 <- bind_rows(lapply(names(esp), function(n) {
  m <- esp[[n]]; ct <- coeftable(m); k <- grep("riego", rownames(ct))
  wr2 <- tryCatch(as.numeric(r2(m, "wr2")), error = function(e) NA_real_)
  tibble(especificacion = n, termino = rownames(ct)[k], b = ct[k, 1], ee = ct[k, 2], p = ct[k, 4],
         efecto_pct = 100 * (exp(ct[k, 1]) - 1), n = nobs(m), r2_within = wr2)
}))
write_csv(t9, file.path(dir_tab, "T9_robustez_riego.csv"))
print(t9, n = Inf, width = Inf)

# ---- Figura 5: coeficientes por especificación + efecto marginal según déficit --------------
g1 <- t9 |> filter(termino %in% c("riego", "riego_tec", "riego_grav")) |>
  mutate(etq = ifelse(termino == "riego", especificacion, paste0(especificacion, ": ", termino)),
         etq = factor(etq, levels = rev(unique(etq)))) |>
  ggplot(aes(efecto_pct, etq)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = 100 * (exp(b - 1.96 * ee) - 1), xmax = 100 * (exp(b + 1.96 * ee) - 1)), height = .2) +
  geom_point(size = 2, colour = "#b2182b") +
  labs(x = "Efecto del riego sobre el rendimiento (%; IC 95 %)", y = NULL, title = "A. Riego según especificación") +
  theme_minimal(base_size = 9)
m5 <- esp[[5]]; b <- coef(m5); V <- vcov(m5)
zz <- seq(quantile(p$z_def, .05, na.rm = TRUE), quantile(p$z_def, .95, na.rm = TRUE), length.out = 50)
me <- tibble(def_mm = zz * sd(p$def_d, na.rm = TRUE) + mean(p$def_d, na.rm = TRUE),
             b = b["riego"] + b["riego:z_def"] * zz,
             se = sqrt(V["riego", "riego"] + zz^2 * V["riego:z_def", "riego:z_def"] + 2 * zz * V["riego", "riego:z_def"]))
g2 <- ggplot(me, aes(def_mm, 100 * (exp(b) - 1))) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = 100 * (exp(b - 1.96 * se) - 1), ymax = 100 * (exp(b + 1.96 * se) - 1)), fill = "#2166ac", alpha = .2) +
  geom_line(colour = "#2166ac", linewidth = 1) +
  labs(x = "Déficit hídrico de la campaña (mm)", y = "Efecto del riego (%)",
       title = "B. Efecto marginal del riego según déficit (esp. 5)") + theme_minimal(base_size = 9)
ggsave(file.path(dir_fig, "F5_robustez_riego.png"), g1 + g2 + plot_layout(widths = c(1.3, 1)),
       width = 11, height = 4.5, dpi = 300, bg = "white")
message("Robustez del riego lista.")
