# =============================================================================
# 10_mgwr_inferencia.R — respuesta a revisores (comentario 3)
# Reestima GWR y MGWR con el paquete de referencia mgwr (Python; Oshan et al. 2019),
# que entrega la inferencia exacta: ENP por covariable (traza de las matrices hat
# parciales), t críticos ajustados (da Silva & Fotheringham 2016), intervalos de
# confianza de los anchos de banda y número de condición local (Comber et al. 2023).
# Se compara con la estimación de GWmodel (script 06).
# Requiere: Python con 'pip install mgwr'.
# Salidas: output/tablas/T16_MGWR_inferencia.csv, T16b_MGWR_diagnostico.csv,
#          output/figuras/F4b_MGWR_coeficientes.png, F9_condicion_local.png
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(ggplot2) })

md  <- readRDS(file.path(dir_interim, "modelos_distritales.rds"))
dsf <- st_read(file.path(dir_interim, "distritos_papa_lisa.gpkg"), quiet = TRUE) |>
  mutate(ln_tviaje = log1p(tviaje_d))
Wl  <- readRDS(file.path(dir_interim, "W_distritos.rds"))
X <- md$X; etq <- c(Intercept = "Intercepto", md$etiquetas)
d <- st_drop_geometry(dsf); d[X] <- scale(d[X])
ent <- tibble(ubigeo = d$ubigeo, y = d$y, X_km = Wl$coords[, 1] / 1000, Y_km = Wl$coords[, 2] / 1000) |>
  bind_cols(as_tibble(d[X]))
dir_py <- file.path(dir_interim, "mgwr_py"); dir.create(dir_py, showWarnings = FALSE)
write_csv(ent, file.path(dir_py, "entrada.csv"))

py <- Sys.which("python"); if (!nzchar(py)) py <- Sys.which("python3")
st <- system2(py, c(shQuote("py/mgwr_inferencia.py"), shQuote(file.path(dir_py, "entrada.csv")),
                    shQuote(dir_py)), stdout = TRUE, stderr = TRUE)
cat(st, sep = "\n")

glob <- read_csv(file.path(dir_py, "mgwr_global.csv"), show_col_types = FALSE)
loc  <- read_csv(file.path(dir_py, "mgwr_local.csv"), show_col_types = FALSE, col_types = cols(ubigeo = "c"))
diag <- read_csv(file.path(dir_py, "mgwr_diag.csv"), show_col_types = FALSE)
stopifnot(identical(as.integer(loc$ubigeo), as.integer(d$ubigeo))); loc$ubigeo <- d$ubigeo

t16 <- glob |> mutate(etiqueta = etq[variable], bw_GWmodel = md$bws[variable],
  media = sapply(variable, function(v) mean(loc[[v]])), p25 = sapply(variable, function(v) quantile(loc[[v]], .25)),
  mediana = sapply(variable, function(v) median(loc[[v]])), p75 = sapply(variable, function(v) quantile(loc[[v]], .75)),
  pct_sig_pos = sapply(seq_along(variable), function(i) 100 * mean(loc[[paste0(variable[i], "_TV")]] >  t_critico[i])),
  pct_sig_neg = sapply(seq_along(variable), function(i) 100 * mean(loc[[paste0(variable[i], "_TV")]] < -t_critico[i]))) |>
  relocate(etiqueta, variable)
write_csv(t16, file.path(dir_tab, "T16_MGWR_inferencia.csv"))
write_csv(diag, file.path(dir_tab, "T16b_MGWR_diagnostico.csv"))
print(t16, width = Inf); print(diag)

# Moran de los residuos MGWR (mgwr)
suppressPackageStartupMessages(library(spdep))
mi <- moran.mc(loc$residual, Wl$W$knn6, nsim = 999)
message("Moran residuos MGWR (mgwr): I = ", round(mi$statistic, 3), ", p = ", mi$p.value)
write_csv(tibble(I = unname(mi$statistic), p = mi$p.value), file.path(dir_tab, "T16c_MGWR_moran_residuos.csv"))

# Figura 4b: coeficientes con t crítico exacto
peru <- leer_distritos() |> st_transform(CRS_M) |> st_union()
dm <- st_transform(dsf, CRS_M)
var_fig <- c("riego", "fert", "semilla_cert", "credito", "def_d", "tmin_d", "meses_helada_d", "ln_tviaje")
tc <- setNames(glob$t_critico, glob$variable); bw <- setNames(glob$bw, glob$variable)
largo <- bind_rows(lapply(var_fig, function(v) st_sf(
  panel = sprintf("%s (bw = %d)", etq[v], round(bw[v])),
  b = ifelse(abs(loc[[paste0(v, "_TV")]]) > tc[v], loc[[v]], NA), geometry = st_geometry(dm))))
largo$panel <- factor(largo$panel, levels = unique(largo$panel))
lim <- max(abs(largo$b), na.rm = TRUE)
f4 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey70", linewidth = .15) +
  geom_sf(data = largo, aes(fill = b), colour = NA) + facet_wrap(~ panel, ncol = 4) +
  scale_fill_gradient2("β local\n(z)", low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
                       limits = c(-lim, lim), na.value = "grey85") +
  labs(caption = "Gris: no significativo con el t crítico ajustado por ENP_j (da Silva & Fotheringham, 2016). bw = vecinos. Estimación: mgwr.") +
  theme_void(base_size = 9) + theme(strip.text = element_text(face = "bold", size = 8, margin = margin(b = 3)))
ggsave(file.path(dir_fig, "F4b_MGWR_coeficientes.png"), f4, width = 11, height = 7, dpi = 300, bg = "white")

# Figura 9: número de condición local
f9 <- ggplot() + geom_sf(data = peru, fill = "white", colour = "grey70", linewidth = .15) +
  geom_sf(data = st_sf(cn = loc$cond_local, geometry = st_geometry(dm)), aes(fill = cn), colour = NA) +
  scale_fill_viridis_c("Número de\ncondición local", option = "C") +
  labs(caption = "Umbral habitual de multicolinealidad local problemática: 30 (Comber et al., 2023).") +
  theme_void(base_size = 9)
ggsave(file.path(dir_fig, "F9_condicion_local.png"), f9, width = 5.5, height = 7, dpi = 300, bg = "white")
message("10 listo.")
