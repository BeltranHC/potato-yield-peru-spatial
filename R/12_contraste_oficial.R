# =============================================================================
# 12_contraste_oficial.R — contraste de la ENA expandida con las estadísticas oficiales
# Cifras oficiales 2024 (MIDAGRI, vía Agencia Andina, 30/05/2025):
#   producción 6,58 Mt; superficie cosechada 337 665 ha; rendimiento 19,5 t/ha.
# Salida: output/tablas/T0c_contraste_oficial.csv
# =============================================================================
source("R/00_config.R")
p <- readRDS(file.path(dir_interim, "ena_papa_ua_todas.rds")) |> filter(ANIO >= 2015)
w <- function(x) ifelse(is.na(x), 0, x)
t0c <- p |> group_by(ANIO) |>
  summarise(n_ua = n(),
            rend_mediano_UA_t = median(rend) / 1000,
            rend_medio_UA_pond_t = weighted.mean(rend, w(FACTOR)) / 1000,
            produccion_expandida_Mt = sum(prod_kg * w(FACTOR)) / 1e9,
            superficie_expandida_ha = sum(sup_cos * w(FACTOR)),
            .groups = "drop") |>
  mutate(rend_agregado_t = produccion_expandida_Mt * 1e6 / superficie_expandida_ha,
         oficial_produccion_Mt = ifelse(ANIO == 2024, 6.58, NA),
         oficial_superficie_ha = ifelse(ANIO == 2024, 337665, NA),
         oficial_rend_t = ifelse(ANIO == 2024, 19.5, NA))
write_csv(t0c, file.path(dir_tab, "T0c_contraste_oficial.csv"))
print(t0c, width = Inf)
