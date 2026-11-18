# =============================================================================
# 01_descargar_ENA2025.R — descarga y descomprime los 21 módulos ENA 2025 (SPSS)
# Fuente: INEI, Microdatos > ENA > 2025 > Anual (encuesta 1036), módulos 1893–1913
# Es re-ejecutable: no vuelve a bajar los zip que ya existen y están completos.
# =============================================================================
source("R/00_config.R")

descargar_modulo <- function(m) {
  url  <- sprintf(url_ena25, m)
  zipf <- file.path(dir_zip25, sprintf("1036-Modulo%d.zip", m))
  ok_zip <- function(f) file.exists(f) && !inherits(try(unzip(f, list = TRUE), silent = TRUE), "try-error")
  if (!ok_zip(zipf)) {
    message("Descargando módulo ", m, " ...")
    for (intento in 1:3) {
      r <- try(download.file(url, zipf, mode = "wb", quiet = TRUE), silent = TRUE)
      if (!inherits(r, "try-error") && ok_zip(zipf)) break
      Sys.sleep(5 * intento)
    }
    if (!ok_zip(zipf)) { warning("Falló la descarga del módulo ", m); return(NULL) }
  } else message("Módulo ", m, ": zip ya presente")

  # Descomprimir en su propia subcarpeta (data/raw/ENA_2025/sav/<modulo>/)
  dest <- file.path(dir_sav25, m)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  cont <- unzip(zipf, list = TRUE)
  unzip(zipf, exdir = dest, junkpaths = TRUE)
  tibble(modulo = m, archivo = basename(cont$Name), bytes = cont$Length)
}

inventario <- bind_rows(lapply(modulos25, descargar_modulo)) |>
  filter(bytes > 0)

write_csv(inventario, file.path(dir_raw25, "inventario_modulos_2025.csv"))
print(inventario |> group_by(modulo) |>
        summarise(n_archivos = n(),
                  sav = paste(archivo[grepl("\\.sav$", archivo, TRUE)], collapse = "; "),
                  MB = round(sum(bytes) / 2^20, 1)), n = Inf)
message("Listo. Revisar también los PDF/diccionarios dentro de ", dir_sav25)
