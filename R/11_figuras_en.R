# =============================================================================
# 11_figuras_en.R — versiones en inglés de las figuras del manuscrito
# Rehace F1–F9 desde los resultados guardados (no reestima los modelos principales;
# solo reajusta modelos rápidos de fixest para los errores estándar de F5 y F8 y
# recalcula las distancias de CAST::geodist para F6).
# Salida: output/figuras_en/*.png (mismos nombres que en output/figuras/)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(ggplot2); library(patchwork); library(fixest)
  library(gstat); library(ranger); library(CAST) })
set.seed(2025)
dir_fig_en <- file.path(dir_out, "figuras_en"); dir.create(dir_fig_en, showWarnings = FALSE)
guardar <- function(p, f, w, h) ggsave(file.path(dir_fig_en, f), p, width = w, height = h, dpi = 300, bg = "white")

var_en <- c(Intercept = "Intercept", riego = "Irrigation", semilla_cert = "Certified seed",
            fert = "Chemical fertiliser", asist = "Technical assistance", credito = "Credit",
            asoc = "Association membership", ln_sup = "ln harvested area", def_d = "Water deficit",
            tmin_d = "Tmin", meses_helada_d = "Frost months", pend_d = "Slope", ln_tviaje = "ln travel time to city")
peru <- leer_distritos() |> st_transform(CRS_M) |> st_union()
dsf  <- st_read(file.path(dir_interim, "distritos_papa_lisa.gpkg"), quiet = TRUE) |> st_transform(CRS_M)
tema <- theme_void(base_size = 10) + theme(legend.position = c(0.2, 0.25), plot.title = element_text(face = "bold"))

# ---- F1: yield map ---------------------------------------------------------------------------
dsf$rend_geo <- exp(dsf$y) / 1000
f1 <- ggplot() + geom_sf(data = peru, fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = dsf, aes(fill = cut(rend_geo, quantile(rend_geo, 0:5 / 5), include.lowest = TRUE)), colour = NA) +
  scale_fill_viridis_d("t/ha (quintiles)") +
  labs(title = "Potato yield by district, 2015–2025",
       caption = sprintf("Expansion-factor weighted geometric mean; %d districts with ≥ 15 FUs. Source: ENA-INEI.", nrow(dsf))) + tema
guardar(f1, "F1_rendimiento_distrital.png", 6, 7.5)

# ---- F2: LISA ---------------------------------------------------------------------------------
lis_en <- c("Alto-Alto" = "High-High", "Bajo-Bajo" = "Low-Low", "Alto-Bajo" = "High-Low",
            "Bajo-Alto" = "Low-High", "No significativo" = "Not significant")
dsf$cl_en <- factor(lis_en[as.character(dsf$cluster)], levels = unname(lis_en))
f2 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = dsf, aes(fill = cl_en), colour = NA) +
  scale_fill_manual("LISA (FDR q < 0.05)", values = c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
    "High-Low" = "#fdae61", "Low-High" = "#abd9e9", "Not significant" = "grey85"), drop = FALSE) +
  labs(title = "Local clusters of potato log yield",
       caption = "Local Moran's I, 9,999 permutations, W = kNN-6; Benjamini–Hochberg correction.") + tema
guardar(f2, "F2_LISA_FDR.png", 6, 7.5)

# ---- F3: variogram ----------------------------------------------------------------------------
vgi <- readRDS(file.path(dir_interim, "variograma_ua.rds"))
lin <- function(vm, tag) { d <- variogramLine(vm, maxdist = 300); d$serie <- tag; d }
s1 <- "Log yield"; s2 <- "OLS residual (environment)"
pts <- bind_rows(mutate(as.data.frame(vgi$vg_y)[c("dist", "gamma", "np")], serie = s1),
                 mutate(as.data.frame(vgi$vg_res)[c("dist", "gamma", "np")], serie = s2))
f3 <- ggplot(pts, aes(dist, gamma, colour = serie)) + geom_point(aes(size = np), alpha = .7) +
  geom_line(data = bind_rows(lin(vgi$vm_y, s1), lin(vgi$vm_res, s2))) +
  scale_size_continuous("Pairs", range = c(0.8, 3)) + scale_colour_manual(NULL, values = c("#1b9e77", "#d95f02")) +
  labs(x = "Distance (km)", y = "Semivariance", title = "Empirical and fitted variograms (FUs 2023–2024)",
       caption = sprintf("Segment-level coordinates. Practical range of the residual ≈ %.0f km.", vgi$rango_bloque_km)) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
guardar(f3, "F3_variograma_UA.png", 7, 4.5)

# ---- F4b and F9: MGWR (mgwr) ---------------------------------------------------------------------
dir_py <- file.path(dir_interim, "mgwr_py")
glob <- read_csv(file.path(dir_py, "mgwr_global.csv"), show_col_types = FALSE)
loc  <- read_csv(file.path(dir_py, "mgwr_local.csv"), show_col_types = FALSE)
tc <- setNames(glob$t_critico, glob$variable); bw <- setNames(glob$bw, glob$variable)
var_fig <- c("riego", "fert", "semilla_cert", "credito", "def_d", "tmin_d", "meses_helada_d", "ln_tviaje")
largo <- bind_rows(lapply(var_fig, function(v) st_sf(panel = sprintf("%s (bw = %d)", var_en[v], round(bw[v])),
  b = ifelse(abs(loc[[paste0(v, "_TV")]]) > tc[v], loc[[v]], NA), geometry = st_geometry(dsf))))
largo$panel <- factor(largo$panel, levels = unique(largo$panel))
lim <- max(abs(largo$b), na.rm = TRUE)
f4 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey70", linewidth = .15) +
  geom_sf(data = largo, aes(fill = b), colour = NA) + facet_wrap(~ panel, ncol = 4) +
  scale_fill_gradient2("Local β\n(z)", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey85") +
  labs(caption = "Grey: not significant at the critical t adjusted by ENP_j (da Silva & Fotheringham, 2016). bw = neighbours. Estimation: mgwr.") +
  theme_void(base_size = 9) + theme(strip.text = element_text(face = "bold", size = 8, margin = margin(b = 3)))
guardar(f4, "F4b_MGWR_coeficientes.png", 11, 7)
f9 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey70", linewidth = .15) +
  geom_sf(data = st_sf(cn = loc$cond_local, geometry = st_geometry(dsf)), aes(fill = cn), colour = NA) +
  scale_fill_viridis_c("Local\ncondition\nnumber", option = "C") +
  labs(caption = "Usual threshold for problematic local multicollinearity: 30 (Comber et al., 2023).") + theme_void(base_size = 9)
guardar(f9, "F9_condicion_local.png", 5.5, 7)

# ---- F5 and F8: FU-level fixed-effects models ------------------------------------------------------
pct <- function(b) 100 * (exp(b) - 1)
u <- readRDS(file.path(dir_interim, "papa_ua_cov.rds")) |> filter(ANIO >= 2015, !is.na(ubigeo_mapa)) |>
  mutate(mes_siembra = ifelse(ventana_declarada, as.integer(format(ini, "%m")), 0L),
         cred_na = as.integer(is.na(credito)),
         across(c(credito, semilla_cert, fert, abono, plag, asist, asoc), ~ ifelse(is.na(.x), 0, .x)),
         z_def = as.numeric(scale(def_d)), z_tmin = as.numeric(scale(tmin_d)), z_hel = as.numeric(scale(meses_helada_d)),
         segmento = sub("^([0-9]+_[^_]+)_.*$", "\\1", id_ua))
ctrl <- "semilla_cert + fert + abono + plag + asist + credito + cred_na + asoc + ln_sup"
m5 <- feols(as.formula(paste("y ~ riego + riego:z_def + riego:z_hel +", ctrl,
                             "+ z_def + z_tmin + z_hel | ubigeo_mapa + ANIO + mes_siembra")), u, cluster = ~ubigeo_mapa)
t9 <- read_csv(file.path(dir_tab, "T9_robustez_riego.csv"), show_col_types = FALSE)
esp_en <- c("1" = "1 Pooled, no controls", "2" = "2 + management + year FE", "3" = "3 + district FE (within district)",
            "4" = "4 + season climate + sowing-month FE", "5" = "5 + irrigation × deficit, irrigation × frost",
            "6" = "6 Technified vs. gravity", "7" = "7 Irrigation × highlands",
            "8" = "8 Weighted by expansion factor", "9" = "9 Segment × year FE (2023–2025, ~1 km)")
g1 <- t9 |> filter(termino %in% c("riego", "riego_tec", "riego_grav")) |>
  mutate(etq = esp_en[substr(especificacion, 1, 1)],
         etq = ifelse(termino == "riego", etq, paste0(etq, ": ", ifelse(termino == "riego_tec", "technified", "gravity"))),
         etq = factor(etq, levels = rev(unique(etq)))) |>
  ggplot(aes(efecto_pct, etq)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = pct(b - 1.96 * ee), xmax = pct(b + 1.96 * ee)), width = .2, orientation = "y") +
  geom_point(size = 2, colour = "#b2182b") +
  labs(x = "Effect of irrigation on yield (%; 95% CI)", y = NULL, title = "A. Irrigation by specification") +
  theme_minimal(base_size = 9)
b <- coef(m5); V <- vcov(m5)
zz <- seq(quantile(u$z_def, .05, na.rm = TRUE), quantile(u$z_def, .95, na.rm = TRUE), length.out = 50)
me <- tibble(def_mm = zz * sd(u$def_d, na.rm = TRUE) + mean(u$def_d, na.rm = TRUE),
             b = b["riego"] + b["riego:z_def"] * zz,
             se = sqrt(V["riego", "riego"] + zz^2 * V["riego:z_def", "riego:z_def"] + 2 * zz * V["riego", "riego:z_def"]))
g2 <- ggplot(me, aes(def_mm, pct(b))) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_ribbon(aes(ymin = pct(b - 1.96 * se), ymax = pct(b + 1.96 * se)), fill = "#2166ac", alpha = .2) +
  geom_line(colour = "#2166ac", linewidth = 1) +
  labs(x = "Seasonal climatic water deficit (mm)", y = "Effect of irrigation (%)",
       title = "B. Marginal effect of irrigation by deficit (spec. 5)") + theme_minimal(base_size = 9)
guardar(g1 + g2 + plot_layout(widths = c(1.3, 1)), "F5_robustez_riego.png", 11, 4.5)

manejo <- c("riego", "semilla_cert", "fert", "asist", "credito", "asoc", "ln_sup")
rhs <- "riego + semilla_cert + fert + asist + credito + asoc + ln_sup + abono + plag + cred_na"
m_ua  <- feols(as.formula(paste("y ~", rhs, "+ z_def + z_tmin + z_hel | ubigeo_mapa + ANIO + mes_siembra")), u, cluster = ~ubigeo_mapa)
m_seg <- feols(as.formula(paste("y ~", rhs, "| segmento + mes_siembra")), filter(u, ANIO >= 2023), cluster = ~ubigeo_mapa)
t11 <- read_csv(file.path(dir_tab, "T11_consistencia_escalas.csv"), show_col_types = FALSE)
imp <- read_csv(file.path(dir_tab, "T7_impactos_SDM_SDEM.csv"), show_col_types = FALSE) |> filter(modelo == "SDM", efecto == "direct")
dz <- st_drop_geometry(dsf) |> mutate(ln_tviaje = log1p(tviaje_d))
sd_x <- sapply(manejo, function(v) sd(dz[[v]], na.rm = TRUE))
se_sdm <- imp$estimacion[match(manejo, imp$variable)] / imp$z[match(manejo, imp$variable)] / sd_x
esc <- c("Between districts (SDM, direct effect)", "Within district (FU, district FE)", "Within segment (FU 2023–2025, ~1 km)")
fd <- bind_rows(
  tibble(etiqueta = var_en[manejo], escala = esc[1], b = t11$sdm_b[match(manejo, t11$variable)], ee = se_sdm),
  tibble(etiqueta = var_en[manejo], escala = esc[2], b = coef(m_ua)[manejo], ee = se(m_ua)[manejo]),
  tibble(etiqueta = var_en[manejo], escala = esc[3], b = coef(m_seg)[manejo], ee = se(m_seg)[manejo]))
fd$etiqueta <- factor(fd$etiqueta, levels = rev(var_en[manejo])); fd$escala <- factor(fd$escala, levels = esc)
f8 <- ggplot(fd, aes(pct(b), etiqueta, colour = escala)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbar(aes(xmin = pct(b - 1.96 * ee), xmax = pct(b + 1.96 * ee)), width = .25,
                position = position_dodge(width = .6), orientation = "y") +
  geom_point(size = 2, position = position_dodge(width = .6)) +
  scale_colour_manual(NULL, values = c("#b2182b", "#2166ac", "#1b7837")) +
  labs(x = "Effect on yield (%; adoption 0 → 1; 95% CI)", y = NULL, title = "Cross-scale consistency of management effects") +
  coord_cartesian(xlim = c(-60, 120)) + theme_minimal(base_size = 9) + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))
guardar(f8, "F8_consistencia_escalas.png", 7.5, 5)

# ---- F6: nearest-neighbour distances ---------------------------------------------------------------
source("R/funciones_prediccion.R")
pt <- preparar_puntos(readRDS(file.path(dir_interim, "papa_ua_cov.rds")))
tr <- pt[pt$ANIO %in% 2023:2024, ]; te <- pt[pt$ANIO == 2025, ]
folds <- readRDS(file.path(dir_interim, "folds_ua_2023_2024.rds"))
dominio <- st_read(file.path(dir_interim, "distritos_papa.gpkg"), quiet = TRUE) |>
  st_transform(CRS_M) |> st_union() |> st_make_valid() |> st_simplify(dTolerance = 500)
cv_en <- c(aleatoria = "random", bloques = "spatial blocks", kNNDM = "kNNDM")
gd <- bind_rows(lapply(names(cv_en), function(esq) {
  g <- tryCatch(geodist(tr, modeldomain = dominio, CVtest = folds[[esq]], testdata = te, samplesize = 2000),
                error = function(e) NULL)
  if (is.null(g)) return(NULL)
  as.data.frame(g) |> mutate(esquema = cv_en[esq])
})) |>
  mutate(what = as.character(what),
         serie = case_when(what == "sample-to-sample" ~ "sample-to-sample",
                           what == "prediction-to-sample" ~ "prediction domain-to-sample",
                           what == "test-to-sample" ~ "2025 test-to-sample",
                           grepl("CV", what) ~ paste0("CV ", esquema), TRUE ~ what)) |>
  distinct(serie, dist, .keep_all = TRUE)
f6 <- ggplot(gd, aes(dist / 1000, colour = serie)) + geom_density(linewidth = .8) +
  scale_x_continuous("Nearest-neighbour distance to training data (km)", trans = "log1p",
                     breaks = c(0, 1, 5, 10, 25, 50, 100, 250)) +
  labs(y = "Density", colour = NULL, title = "Does cross-validation mimic the actual prediction distance?") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
guardar(f6, "F6_distancias_CV.png", 8, 4.8)

# ---- F7: prediction map with AOA -----------------------------------------------------------------------
mapa <- readRDS(file.path(dir_interim, "mapa_RK_2025.rds"))
mapa_in <- mapa[mapa$AOA, ]; mapa_out <- mapa[!mapa$AOA, ]
base <- function() list(geom_sf(data = peru, fill = "white", colour = "grey60", linewidth = .2), theme_void(base_size = 9))
fuera <- geom_sf(data = mapa_out, colour = "grey80", size = .25, shape = 15)
f7 <- (ggplot() + base() + fuera + geom_sf(data = mapa_in, aes(colour = rend_t), size = .25, shape = 15) +
         scale_colour_viridis_c("t/ha") + labs(title = "A. Predicted yield (RK, environment)")) +
      (ggplot() + base() + fuera + geom_sf(data = mapa_in, aes(colour = sd_krig), size = .25, shape = 15) +
         scale_colour_viridis_c("Kriging SD\n(log)", option = "magma") + labs(title = "B. Kriging uncertainty")) +
      plot_annotation(caption = paste0("5 km grid over potato-growing districts; October 2024–April 2025 season. ",
                                       "Grey: outside the area of applicability (AOA; Meyer & Pebesma, 2021)."))
guardar(f7, "F7_mapa_RK_2025_AOA.png", 10, 6.5)
message("11 listo: ", length(list.files(dir_fig_en, "[.]png$")), " figuras en ", dir_fig_en)
