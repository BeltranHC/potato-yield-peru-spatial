# =============================================================================
# 04_covariables.R — covariables ambientales y de acceso para UA y distritos
# Fuentes públicas (no requieren cuenta):
#   * Límites distritales con UBIGEO (IDDIST): peru-geojson (base INEI)
#   * Relieve: SRTM 90 m (CGIAR, vía geodata) -> altitud y pendiente
#   * Acceso: tiempo de viaje a ciudad >= 20 mil hab. (Nelson et al. 2019) — reemplaza
#     a P225 (horas a la capital distrital), que la ENA solo pregunta hasta 2019
#   * Clima mensual: TerraClimate (~4 km, 1958–2025): ppt, tmin, tmax, def (déficit hídrico)
# Opcional (si se colocan los archivos): PISCOp v2.1 (IRI, requiere cuenta) para robustez.
#
# Ventana climática de cada UA = meses entre siembra (P209) y cosecha (P207 fin)
# declarados. Si faltan o son inválidas (> 12 meses), se usa la campaña andina típica:
# octubre del año anterior a abril del año de la ENA.
#
# Salidas:
#   data/interim/papa_ua_cov.rds   UA-campaña con covariables (clima distrital para todas;
#                                  covariables puntuales *_pt para UA con coordenadas)
#   data/interim/distritos_papa.gpkg  polígonos + agregados distritales (para 05 y 06)
# =============================================================================
source("R/00_config.R")
suppressPackageStartupMessages({ library(sf); library(terra); library(exactextractr) })
instalar_faltantes("geodata")

BUF_PT   <- 1000   # m: las coordenadas 2023–2025 son del segmento (2025 redondeadas a 0,01°)
ANIOS_CL <- 2013:2025
VARS_CL  <- c("ppt", "tmin", "tmax", "def")
T_BASE   <- 4      # °C, temperatura base para grados-día de papa

dir_lim <- file.path(dir_ext, "limites");  dir_dem <- file.path(dir_ext, "dem")
dir_acc <- file.path(dir_ext, "acceso");   dir_cl  <- file.path(dir_ext, "terraclimate")
for (d in c(dir_lim, dir_dem, dir_acc, dir_cl)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

papa <- readRDS(file.path(dir_interim, "ena_papa_ua_todas.rds"))

# ---- 1. Límites distritales -------------------------------------------------------
f_dist <- file.path(dir_lim, "distritos.gpkg")
if (!file.exists(f_dist)) {
  url <- "https://raw.githubusercontent.com/juaneladio/peru-geojson/master/peru_distrital_simple.geojson"
  st_read(url, quiet = TRUE) |>
    transmute(ubigeo = IDDIST, distrito = NOMBDIST, provincia = NOMBPROV, departamento = NOMBDEP) |>
    st_make_valid() |> st_write(f_dist, quiet = TRUE)
}
dist <- leer_distritos(f_dist)
# Distritos creados después del mapa: se asignan al distrito capital de su provincia (PP01)
sin_poly <- setdiff(unique(papa$ubigeo), dist$ubigeo)
remap <- setNames(ifelse(paste0(substr(sin_poly, 1, 4), "01") %in% dist$ubigeo,
                         paste0(substr(sin_poly, 1, 4), "01"), NA), sin_poly)
papa$ubigeo_mapa <- ifelse(papa$ubigeo %in% names(remap), remap[papa$ubigeo], papa$ubigeo)
message("UBIGEO sin polígono: ", length(sin_poly), " (", sum(papa$ubigeo %in% sin_poly),
        " UA) -> reasignados a la capital provincial")

# ---- 2. Relieve: SRTM 90 m -------------------------------------------------------------
f_dem <- file.path(dir_dem, "srtm90_andes.tif"); f_pend <- file.path(dir_dem, "pendiente90_andes.tif")
if (!file.exists(f_dem)) {
  # teselas de 5° que contienen distritos con papa
  bb <- st_bbox(dist[dist$ubigeo %in% papa$ubigeo_mapa, ])
  lon <- seq(floor(bb["xmin"] / 5) * 5 + 2.5, ceiling(bb["xmax"] / 5) * 5 - 2.5, by = 5)
  lat <- seq(floor(bb["ymin"] / 5) * 5 + 2.5, ceiling(bb["ymax"] / 5) * 5 - 2.5, by = 5)
  tiles <- list()
  for (x in lon) for (y in lat) {
    r <- try(geodata::elevation_3s(lon = x, lat = y, path = dir_dem), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r)) tiles[[length(tiles) + 1]] <- r
  }
  dem <- crop(merge(sprc(tiles)), ext(bb) + 0.1, filename = f_dem, overwrite = TRUE)
  terrain(dem, "slope", unit = "degrees", filename = f_pend, overwrite = TRUE)
}
topo <- c(rast(f_dem), rast(f_pend)); names(topo) <- c("alt", "pend")

# ---- 3. Acceso a mercado -----------------------------------------------------------------
acc <- geodata::travel_time(to = "city", size = 7, up = TRUE, path = dir_acc)
acc <- crop(acc, ext(topo)); names(acc) <- "tviaje_min"

# ---- 4. Clima mensual TerraClimate (recorte Perú vía NCSS) -------------------------------
url_tc <- paste0("http://thredds.northwestknowledge.net:8080/thredds/ncss/TERRACLIMATE_ALL/data/",
                 "TerraClimate_%s_%d.nc?var=%s&north=0.5&west=-81.5&east=-68.5&south=-18.5",
                 "&disableProjSubset=on&horizStride=1&temporal=all&accept=netcdf")
for (v in VARS_CL) for (a in ANIOS_CL) {
  f <- file.path(dir_cl, sprintf("tc_%s_%d.nc", v, a))
  if (!file.exists(f) || file.size(f) < 1e6) {
    message("TerraClimate ", v, " ", a)
    for (intento in 1:5) {
      r <- try(download.file(sprintf(url_tc, v, a, v), f, mode = "wb", quiet = TRUE), silent = TRUE)
      if (!inherits(r, "try-error") && file.size(f) > 1e6) break
      message("  reintento ", intento); Sys.sleep(10 * intento)
    }
    if (!file.exists(f) || file.size(f) < 1e6) stop("No se pudo descargar TerraClimate ", v, " ", a)
  }
}
clima <- lapply(setNames(VARS_CL, VARS_CL), function(v) {
  r <- rast(file.path(dir_cl, sprintf("tc_%s_%d.nc", v, ANIOS_CL)))
  time(r) <- seq(as.Date(sprintf("%d-01-01", min(ANIOS_CL))), by = "month", length.out = nlyr(r))
  names(r) <- paste0(v, "_", format(time(r), "%Y%m"))
  r
})
meses <- as.Date(time(clima$ppt))
dias_mes <- as.numeric(diff(c(meses, seq(max(meses), by = "month", length.out = 2)[2])))

# ---- 5. Ventana fenológica de cada UA ----------------------------------------------------
papa <- papa |>
  mutate(ini = as.Date(ifelse(is.finite(siembra), siembra, NA)),
         fin = as.Date(ifelse(is.finite(cosecha), cosecha, NA)),
         ok  = !is.na(ini) & !is.na(fin) & fin >= ini &
               as.numeric(fin - ini) <= 366 & ini >= min(meses) & fin <= max(meses),
         ventana_declarada = ok,
         ini = if_else(ok, ini, as.Date(sprintf("%d-10-01", ANIO - 1))),
         fin = if_else(ok, fin, as.Date(sprintf("%d-04-01", ANIO)))) |>
  select(-ok)
message("UA con ventana declarada válida: ", round(100 * mean(papa$ventana_declarada), 1), " %")

# Resume una matriz n x meses en la ventana [ini, fin] de cada fila
resumir_ventana <- function(M, ini, fin) {
  i0 <- match(ini, meses); i1 <- match(fin, meses)
  out <- matrix(NA_real_, length(ini), 6,
                dimnames = list(NULL, c("ppt", "tmin", "tmax", "def", "meses_helada", "gdd")))
  for (k in which(!is.na(i0) & !is.na(i1))) {
    j <- i0[k]:i1[k]
    tmn <- M$tmin[k, j]; tmx <- M$tmax[k, j]
    out[k, ] <- c(sum(M$ppt[k, j]), mean(tmn), mean(tmx), sum(M$def[k, j]),
                  sum(tmn < 0), sum(pmax(0, (tmn + tmx) / 2 - T_BASE) * dias_mes[j]))
  }
  out
}

# ---- 6. Clima a nivel distrital (todas las UA, con o sin coordenadas) --------------------
dist_p <- dist[dist$ubigeo %in% papa$ubigeo_mapa, ]
f_cd <- file.path(dir_interim, "clima_distrito_mes.rds")
if (!file.exists(f_cd)) {
  cd <- lapply(clima, function(r) as.matrix(exact_extract(r, dist_p, "mean", progress = FALSE)))
  saveRDS(list(ubigeo = dist_p$ubigeo, M = cd), f_cd)
}
cd <- readRDS(f_cd)
idx <- match(papa$ubigeo_mapa, cd$ubigeo)
Mua <- lapply(cd$M, function(m) m[idx, , drop = FALSE])
cl_d <- resumir_ventana(Mua, papa$ini, papa$fin)
colnames(cl_d) <- paste0(colnames(cl_d), "_d")
papa <- bind_cols(papa, as_tibble(cl_d))

# Relieve y acceso distritales (zonal)
zon <- cbind(exact_extract(topo, dist_p, "mean", progress = FALSE),
             exact_extract(acc, dist_p, "mean", progress = FALSE))
dist_p[c("alt_d", "pend_d", "tviaje_d")] <- zon

# ---- 7. Covariables puntuales (UA con coordenadas; buffer de 1 km) ----------------------
pt <- papa |> filter(!is.na(LONGITUD), !is.na(LATITUD), between(LATITUD, -19, 1), between(LONGITUD, -82, -68)) |>
  st_as_sf(coords = c("LONGITUD", "LATITUD"), crs = CRS_GEO, remove = FALSE)
buf <- st_buffer(st_transform(pt, CRS_M), BUF_PT) |> st_transform(CRS_GEO)
v_pt <- cbind(exact_extract(topo, buf, "mean", progress = FALSE),
              exact_extract(acc, buf, "mean", progress = FALSE))
names(v_pt) <- c("alt_pt", "pend_pt", "tviaje_pt")
Mpt <- lapply(clima, function(r) as.matrix(exact_extract(r, buf, "mean", progress = FALSE)))
cl_p <- resumir_ventana(Mpt, pt$ini, pt$fin); colnames(cl_p) <- paste0(colnames(cl_p), "_pt")
pt_cov <- bind_cols(st_drop_geometry(pt)["id_ua"], v_pt, as_tibble(cl_p))
papa <- left_join(papa, pt_cov, by = "id_ua") |>
  mutate(qc_alt = ifelse(is.na(alt_pt) | is.na(ALTITUD), NA, abs(alt_pt - ALTITUD) < 500))
message("UA con covariables puntuales: ", nrow(pt_cov))

# Control: ¿el punto cae en el distrito declarado?
en <- st_join(pt["id_ua"], dist["ubigeo"], join = st_within) |> st_drop_geometry() |>
  distinct(id_ua, .keep_all = TRUE) |> rename(ubigeo_pt = ubigeo)
papa <- left_join(papa, en, by = "id_ua") |>
  mutate(qc_dist = ifelse(is.na(ubigeo_pt), NA, ubigeo_pt == ubigeo_mapa))

saveRDS(papa, file.path(dir_interim, "papa_ua_cov.rds"))

# ---- 8. Agregación distrital (pool 2015–2025, pesos = FACTOR) ----------------------------
covs_ua <- c("y", "riego", "riego_tec", "semilla_cert", "fert", "abono", "plag", "asist",
             "credito", "seguro", "asoc", "ln_sup",
             "ppt_d", "tmin_d", "tmax_d", "def_d", "meses_helada_d", "gdd_d")
dist_df <- papa |> filter(ANIO >= 2015, !is.na(ubigeo_mapa)) |>
  mutate(w = ifelse(is.na(FACTOR), 0, FACTOR)) |>
  group_by(ubigeo = ubigeo_mapa) |>
  summarise(n_ua = n(), n_anios = n_distinct(ANIO),
            across(all_of(covs_ua), ~ { ok <- !is.na(.x) & w > 0
                                        if (any(ok)) weighted.mean(.x[ok], w[ok]) else NA_real_ }),
            .groups = "drop")
dsf <- inner_join(dist_p, dist_df, by = "ubigeo")
st_write(dsf, file.path(dir_interim, "distritos_papa.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 9. Resumen --------------------------------------------------------------------------
res <- papa |> group_by(ANIO) |>
  summarise(n = n(), ventana_decl = mean(ventana_declarada), ppt_d = median(ppt_d, na.rm = TRUE),
            tmin_d = median(tmin_d, na.rm = TRUE), helada_d = mean(meses_helada_d > 0, na.rm = TRUE),
            con_pt = mean(!is.na(alt_pt)), qc_dist = mean(qc_dist, na.rm = TRUE),
            qc_alt = mean(qc_alt, na.rm = TRUE))
print(res, n = Inf, width = Inf)
write_csv(res, file.path(dir_tab, "T0b_covariables_por_anio.csv"))
message("Distritos con papa (pool 2015–2025): ", nrow(dsf), " | con >= 15 UA: ", sum(dsf$n_ua >= 15))
