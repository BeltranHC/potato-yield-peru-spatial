# =============================================================================
# 06_modelos_distritales.R — Etapa 2: explicación a nivel distrital (pool 2015–2025)
#   OLS -> diagnósticos (VIF, Moran de residuos, LM robustos) ->
#   SLX / SAR / SEM / SDM / SDEM -> impactos LeSage & Pace (SDM) ->
#   GWR y MGWR (GWmodel), con t críticos ajustados (da Silva & Fotheringham 2016)
# Covariables estandarizadas (z) -> coeficientes comparables entre sí y con MGWR.
# Colinealidad: Tmin ~ altitud ~ grados-día (r ≈ 0,95) y ln ppt ~ déficit (r = -0,88):
#   se usan tmin_d y def_d. VIF máximo del conjunto final < 5.
# Salidas: output/tablas/T4–T7, output/figuras/F4, data/interim/modelos_distritales.rds
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(spdep); library(spatialreg)
                                 library(GWmodel); library(ggplot2) })
set.seed(2025)

dsf <- st_read(file.path(dir_interim, "distritos_papa_lisa.gpkg"), quiet = TRUE) |>
  mutate(ln_tviaje = log1p(tviaje_d))
Wl  <- readRDS(file.path(dir_interim, "W_distritos.rds"))
stopifnot(identical(Wl$ubigeo, dsf$ubigeo))
lw  <- Wl$W$knn6

xs <- c(riego = "Riego", semilla_cert = "Semilla certificada", fert = "Fertilizante químico",
        asist = "Asistencia técnica", credito = "Crédito", asoc = "Asociatividad",
        ln_sup = "ln superficie", def_d = "Déficit hídrico", tmin_d = "Tmin",
        meses_helada_d = "Meses con helada", pend_d = "Pendiente", ln_tviaje = "ln tiempo a ciudad")
X  <- names(xs)
d  <- st_drop_geometry(dsf)
d[X] <- scale(d[X])
fml <- reformulate(X, "y")

# ---- 1. OLS y diagnósticos ------------------------------------------------------------
ols <- lm(fml, data = d)
vif <- sapply(X, function(z) 1 / (1 - summary(lm(reformulate(setdiff(X, z), z), d))$r.squared))
bp  <- { e2 <- residuals(ols)^2; a <- lm(reformulate(X, "e2"), data = cbind(d, e2 = e2))
         s <- nrow(d) * summary(a)$r.squared; c(BP = s, p = pchisq(s, length(X), lower.tail = FALSE)) }
mi_ols <- lm.morantest(ols, lw)
lm_t   <- lm.RStests(ols, lw, test = "all")
num1 <- function(z) if (is.null(z) || !length(z)) NA_real_ else as.numeric(z[[1]])
t5 <- tibble(prueba = names(lm_t),
             estadistico = vapply(lm_t, function(t) num1(t$statistic), 0),
             gl = vapply(lm_t, function(t) num1(t$parameter), 0),
             p = vapply(lm_t, function(t) num1(t$p.value), 0))
t5 <- bind_rows(tibble(prueba = "Moran I residuos OLS", estadistico = num1(mi_ols$estimate),
                       gl = NA_real_, p = num1(mi_ols$p.value)),
                tibble(prueba = "Breusch-Pagan", estadistico = num1(bp["BP"]), gl = length(X), p = num1(bp["p"])),
                t5)
write_csv(t5, file.path(dir_tab, "T5_diagnosticos_OLS_LM.csv"))
print(t5); message("VIF máx = ", round(max(vif), 2))

# ---- 2. Modelos de econometría espacial -----------------------------------------------
slx  <- lmSLX(fml, d, lw)
sar  <- lagsarlm(fml, d, lw)
sem  <- errorsarlm(fml, d, lw)
sdm  <- lagsarlm(fml, d, lw, Durbin = TRUE)
sdem <- errorsarlm(fml, d, lw, Durbin = TRUE)
mods <- list(OLS = ols, SLX = slx, SAR = sar, SEM = sem, SDM = sdm, SDEM = sdem)

res_moran <- function(m) { r <- residuals(m); mm <- moran.mc(r, lw, nsim = 999); c(I = unname(mm$statistic), p = mm$p.value) }
t6 <- tibble(modelo = names(mods),
             logLik = sapply(mods, function(m) as.numeric(logLik(m))),
             k = sapply(mods, function(m) attr(logLik(m), "df")),
             AIC = sapply(mods, AIC),
             rho = c(NA, NA, sar$rho, NA, sdm$rho, NA),
             lambda = c(NA, NA, NA, sem$lambda, NA, sdem$lambda),
             I_residuos = sapply(mods, function(m) res_moran(m)["I"]),
             p_I_residuos = sapply(mods, function(m) res_moran(m)["p"]))
lr <- list(`SDM vs SAR (H0: θ = 0)` = LR.Sarlm(sdm, sar), `SDM vs SEM (H0: θ = -ρβ)` = LR.Sarlm(sdm, sem),
           `SDM vs SLX (H0: ρ = 0)` = LR.Sarlm(sdm, slx), `SDEM vs SEM (H0: θ = 0)` = LR.Sarlm(sdem, sem))
t6_lr <- tibble(contraste = names(lr), LR = sapply(lr, function(t) unname(t$statistic)),
                gl = sapply(lr, function(t) unname(t$parameter)), p = sapply(lr, function(t) t$p.value))

# ---- 3. Impactos (LeSage & Pace 2009) ----------------------------------------------------
trMC <- trW(as(lw, "CsparseMatrix"), type = "MC")
imp_tab <- function(m, nombre) {
  im <- summary(impacts(m, tr = trMC, R = 2000), zstats = TRUE, short = TRUE)
  bind_rows(lapply(seq_along(c("direct", "indirect", "total")), function(j) {
    k <- c("direct", "indirect", "total")[j]
    est <- im$res[[k]]; se <- im$semat[, j]; q <- im[[paste0(k, "_sum")]]$quantiles
    tibble(modelo = nombre, efecto = k, variable = names(est), estimacion = unname(est),
           ee = unname(se), z = unname(est / se), p = 2 * pnorm(-abs(unname(est / se))),
           ic95_inf = unname(q[, "2.5%"]), ic95_sup = unname(q[, "97.5%"]))
  }))
}
imp_sdm <- imp_tab(sdm, "SDM")
# SDEM: directo = β, indirecto = θ (lag.X), total = β + θ (EE con la covarianza de ambos)
imp_sdem <- {
  b <- coef(sdem); V <- vcov(sdem)
  bind_rows(lapply(X, function(v) {
    lv <- paste0("lag.", v)
    est <- c(direct = b[[v]], indirect = b[[lv]], total = b[[v]] + b[[lv]])
    se  <- sqrt(c(V[v, v], V[lv, lv], V[v, v] + V[lv, lv] + 2 * V[v, lv]))
    tibble(modelo = "SDEM", efecto = names(est), variable = v, estimacion = unname(est),
           z = unname(est / se), p = 2 * pnorm(-abs(unname(est / se))))
  }))
}
t7 <- bind_rows(imp_sdm, imp_sdem) |> mutate(etiqueta = xs[variable])
write_csv(t7, file.path(dir_tab, "T7_impactos_SDM_SDEM.csv"))

# Coeficientes de todos los modelos (Tabla 4)
coef_tab <- function(m, nombre) {
  s <- if (inherits(m, "lm")) summary(m)$coefficients else summary(m)$Coef
  tibble(modelo = nombre, termino = rownames(s), b = s[, 1], ee = s[, 2], p = s[, 4])
}
t4 <- bind_rows(lapply(names(mods), function(n) coef_tab(mods[[n]], n)))

# ---- 4. GWR y MGWR ---------------------------------------------------------------------------
co_km <- Wl$coords / 1000
dsp <- sp::SpatialPointsDataFrame(co_km, d[c("y", X)])
f_gw <- file.path(dir_interim, "gwr_mgwr.rds")           # caché: MGWR tarda varios minutos
if (file.exists(f_gw) && identical(readRDS(f_gw)$X, X)) {
  gw <- readRDS(f_gw)
} else {
  t0 <- Sys.time()
  bw_gwr <- bw.gwr(fml, dsp, approach = "AICc", kernel = "bisquare", adaptive = TRUE)
  gwr    <- gwr.basic(fml, dsp, bw = bw_gwr, kernel = "bisquare", adaptive = TRUE)
  message("GWR: bw = ", bw_gwr, " vecinos")
  mgwr <- gwr.multiscale(fml, dsp, kernel = "bisquare", adaptive = TRUE, criterion = "dCVR",
                         max.iterations = 2000, threshold = 1e-5, bws0 = rep(bw_gwr, length(X) + 1),
                         verbose = FALSE, predictor.centered = rep(FALSE, length(X)))
  message("MGWR: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  gw <- list(X = X, bw_gwr = bw_gwr, gwr = gwr, mgwr = mgwr)
  saveRDS(gw, f_gw)
}
gwr <- gw$gwr; mgwr <- gw$mgwr; bw_gwr <- gw$bw_gwr
message("Componentes MGWR: SDF = ", paste(head(names(mgwr$SDF), 40), collapse = " "),
        " | diagnóstico = ", paste(names(mgwr$GW.diagnostic), collapse = " "))
bws <- setNames(mgwr$GW.arguments$bws, c("Intercept", X))

sdf <- as.data.frame(mgwr$SDF)
# ENP por covariable = traza de su matriz hat parcial; aproximación n/bw si no está disponible
enp <- if (!is.null(mgwr$GW.diagnostic$enp.j)) mgwr$GW.diagnostic$enp.j else nrow(d) / bws
alpha_j <- 0.05 / enp
t_crit  <- qt(1 - alpha_j / 2, df = nrow(d) - 1)
names(t_crit) <- names(bws)
resumen_mgwr <- bind_rows(lapply(names(bws), function(v) {
  b <- sdf[[v]]; tv <- sdf[[paste0(v, "_TV")]]
  tibble(variable = v, etiqueta = c(Intercept = "Intercepto", xs)[v], bw_vecinos = bws[v],
         bw_pct = 100 * bws[v] / nrow(d), t_critico = t_crit[v],
         media = mean(b), de = sd(b), min = min(b), p25 = quantile(b, .25), mediana = median(b),
         p75 = quantile(b, .75), max = max(b),
         pct_sig_pos = 100 * mean(tv >  t_crit[v]), pct_sig_neg = 100 * mean(tv < -t_crit[v]))
}))
write_csv(resumen_mgwr, file.path(dir_tab, "T8_MGWR_resumen.csv"))
print(resumen_mgwr, width = Inf)

# Comparación global de ajuste (OLS, espaciales, GWR, MGWR)
r2_pseudo <- function(m) cor(d$y, d$y - residuals(m))^2
mgwr_res <- sdf$residual
t6 <- bind_rows(t6 |> mutate(R2 = sapply(mods, r2_pseudo), AICc = NA_real_),
  tibble(modelo = "GWR", AIC = gwr$GW.diagnostic$AIC, AICc = gwr$GW.diagnostic$AICc,
         R2 = gwr$GW.diagnostic$gw.R2,
         I_residuos = unname(moran.mc(gwr$SDF$residual, lw, 999)$statistic),
         p_I_residuos = moran.mc(gwr$SDF$residual, lw, 999)$p.value),
  tibble(modelo = "MGWR", AIC = mgwr$GW.diagnostic$AIC, AICc = mgwr$GW.diagnostic$AICc,
         R2 = mgwr$GW.diagnostic$R2.val,
         I_residuos = unname(moran.mc(mgwr_res, lw, 999)$statistic),
         p_I_residuos = moran.mc(mgwr_res, lw, 999)$p.value))
write_csv(t6, file.path(dir_tab, "T6_comparacion_modelos.csv"))
write_csv(t6_lr, file.path(dir_tab, "T6b_contrastes_LR.csv"))
write_csv(t4, file.path(dir_tab, "T4_coeficientes_modelos.csv"))
print(t6, width = Inf); print(t6_lr)
print(t7 |> filter(modelo == "SDM") |> select(efecto, etiqueta, estimacion, p) |>
        tidyr::pivot_wider(names_from = efecto, values_from = c(estimacion, p)), width = Inf)

saveRDS(list(ols = ols, slx = slx, sar = sar, sem = sem, sdm = sdm, sdem = sdem, gwr = gwr,
             mgwr = mgwr, bws = bws, t_crit = t_crit, X = X, etiquetas = xs, vif = vif),
        file.path(dir_interim, "modelos_distritales.rds"))

# ---- 5. Figura 4: coeficientes MGWR (no significativos en gris) -------------------------------
peru <- leer_distritos() |> st_transform(CRS_M) |> st_union()
dm <- st_transform(dsf, CRS_M)
var_fig <- c("riego", "fert", "semilla_cert", "credito", "def_d", "tmin_d", "meses_helada_d", "ln_tviaje")
largo <- bind_rows(lapply(var_fig, function(v) {
  b <- sdf[[v]]; tv <- sdf[[paste0(v, "_TV")]]
  st_sf(panel = sprintf("%s (bw = %d)", xs[v], round(bws[v])), b = ifelse(abs(tv) > t_crit[v], b, NA),
        geometry = st_geometry(dm))
}))
largo$panel <- factor(largo$panel, levels = unique(largo$panel))
lim <- max(abs(largo$b), na.rm = TRUE)
f4 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = largo, aes(fill = b), colour = NA) + facet_wrap(~ panel, ncol = 4) +
  scale_fill_gradient2("β local\n(z)", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0, limits = c(-lim, lim), na.value = "grey85") +
  labs(title = "MGWR: efectos locales sobre el ln-rendimiento de papa",
       caption = "Gris: no significativo con t crítico ajustado (da Silva & Fotheringham, 2016). bw = vecinos.") +
  theme_void(base_size = 9) + theme(strip.text = element_text(face = "bold", size = 8))
ggsave(file.path(dir_fig, "F4_MGWR_coeficientes.png"), f4, width = 11, height = 7.5, dpi = 300, bg = "white")
message("Modelos distritales listos.")
