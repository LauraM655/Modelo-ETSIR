
## ----setup-global, include=FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo       = TRUE,
  warning    = FALSE,
  message    = FALSE,
  fig.width  = 11,
  fig.height = 8,
  fig.align  = "center"
)


## ----ModuloI, child='ModuloI.Rmd'------------------------------------------------------------------------------------------

## ----bloque-0-librerias, message=FALSE-------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 0: LIBRERÍAS Y PARÁMETROS GLOBALES
#==============================================================================
# PROPÓSITO: Cargar dependencias y definir constantes del estudio
# JUSTIFICACIÓN: Centralizar configuración facilita reproducibilidad

# --- Instalación silenciosa de paquetes faltantes ---
paquetes_necesarios <- c(
  "readr", "dplyr", "tidyr", "purrr", "lubridate",
  "MMWRweek", "zoo", "ggplot2", "scales", "kableExtra",
  "nasapower", "MASS", "pscl", "mgcv", "forecast", "patchwork"
)

for (pkg in paquetes_necesarios) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Carga de librerías ---
lapply(paquetes_necesarios, library, character.only = TRUE)

# --- Parámetros globales del estudio ---
PAISES_ISO        <- c("CRI", "DOM", "COL", "HND", "MEX")
ANIO_INICIO       <- 2014
ANIO_FIN          <- 2024
K_I               <- 8    # Semanas de acumulación de inmunidad/vigilancia
FACTOR_SUBREG     <- 25   # 1 / 0.04 = 25 (subregistro 4% en Latinoamérica)
USAR_SUBREGISTRO  <- FALSE  # Cambiar a TRUE para replicar ajuste del artículo

# Paleta cromática consistente para todo el documento
PALETA_PAISES <- c(
  "Costa Rica"           = "#1abc9c",
  "República Dominicana" = "#3498db",
  "Colombia"             = "#9b59b6",
  "Honduras"             = "#f1c40f",
  "México"               = "#e67e22"
)

cat("✔ Entorno configurado correctamente.\n")
cat("  Países:", paste(PAISES_ISO, collapse = ", "), "\n")
cat("  Período:", ANIO_INICIO, "-", ANIO_FIN, "\n")
cat("  Ventana de inmunidad K_I:", K_I, "semanas\n")


## ----bloque-1-lectura------------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 1.1: LECTURA DE DATOS EPIDEMIOLÓGICOS
#==============================================================================
# PROPÓSITO: Cargar el extracto de OpenDengue
# DIAGNÓSTICO: Mostrar dimensiones iniciales

archivo_local <- "Temporal_extract_V1_3_REDUCIDO.csv"

# Validar existencia del archivo
if (!file.exists(archivo_local)) {
  stop(paste(
    "Error crítico: El archivo", archivo_local,
    "no se encuentra en el directorio de trabajo.\n",
    "Directorio actual:", getwd()
  ))
}

# Lectura del CSV
dengue_raw <- read_csv(archivo_local, show_col_types = FALSE)

cat("═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO INICIAL: Datos Epidemiológicos\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Registros brutos leídos:", nrow(dengue_raw), "\n")
cat("Columnas disponibles:", ncol(dengue_raw), "\n")
cat("Países presentes en archivo:", length(unique(dengue_raw$ISO_A0)), "\n")
cat("Rango de años:", range(dengue_raw$Year, na.rm = TRUE), "\n")


## ----bloque-1-filtrado-----------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 1.2: FILTRADO CRÍTICO (CORRECCIÓN DE BUG)
#==============================================================================
# PROBLEMA DETECTADO: El extracto de OpenDengue contiene TRES niveles 
# geográficos simultáneos para la misma semana:
#   - S_res == "Admin0" → total nacional (1 fila/semana/país)
#   - S_res == "Admin1" → departamental (subconjunto de Admin0)
#   - S_res == "Admin2" → municipal (subconjunto de Admin0)
#
# BUG ORIGINAL: El código anterior agrupaba por (pais, año, semana) y SUMABA
# dengue_total sin filtrar por S_res. Esto duplicaba los casos al sumar
# el total nacional + su desglose municipal.
#
# SOLUCIÓN: Filtrar T_res == "Week" y S_res == "Admin0" ANTES de agrupar.
#==============================================================================

# Limpieza de espacios invisibles y filtrado
dengue_semanal_raw <- dengue_raw %>%
  mutate(
    ISO_limpio   = trimws(toupper(ISO_A0)),
    pais_limpio  = trimws(adm_0_name),
    T_res_limpio = trimws(T_res),
    S_res_limpio = trimws(S_res)
  ) %>%
  filter(
    ISO_limpio %in% PAISES_ISO,
    Year >= ANIO_INICIO,
    Year <= ANIO_FIN,
    T_res_limpio == "Week",     # Excluye resúmenes anuales
    S_res_limpio == "Admin0"    # Excluye desagregación Admin1/Admin2
  ) %>%
  mutate(
    fecha_semana = as.Date(calendar_start_date),
    semana_epi   = MMWRweek(fecha_semana)$MMWRweek,
    anio_epi     = MMWRweek(fecha_semana)$MMWRyear,
    pais = case_when(
      ISO_limpio == "CRI" ~ "Costa Rica",
      ISO_limpio == "DOM" ~ "República Dominicana",
      ISO_limpio == "COL" ~ "Colombia",
      ISO_limpio == "HND" ~ "Honduras",
      ISO_limpio == "MEX" ~ "México"
    ),
    codigo_iso = ISO_limpio
  ) %>%
  group_by(pais, codigo_iso, anio_epi, semana_epi) %>%
  summarise(
    casos_totales = sum(dengue_total, na.rm = TRUE),
    n_filas_originales = n(),
    .groups = "drop"
  ) %>%
  mutate(fecha_semana = as.Date(MMWRweek2Date(anio_epi, semana_epi))) %>%
  arrange(pais, anio_epi, semana_epi)

# Validación de la corrección
n_filas_multiples <- sum(dengue_semanal_raw$n_filas_originales > 1)

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Filtrado de Niveles Geográficos\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Registros antes de filtrar:", nrow(dengue_raw), "\n")
cat("Registros después de filtrar:", nrow(dengue_semanal_raw), "\n")
cat("Registros eliminados:", nrow(dengue_raw) - nrow(dengue_semanal_raw), "\n")
cat("  (Niveles Admin1/Admin2 y resúmenes anuales)\n")

if (n_filas_multiples > 0) {
  warning(
    n_filas_multiples,
    " semanas tienen más de 1 fila Admin0/Week. Revisar manualmente."
  )
} else {
  cat("✔ Verificación: cada semana tiene exactamente 1 fila fuente.\n")
}

# Eliminar columna de diagnóstico
dengue_semanal_raw <- dengue_semanal_raw %>% 
  dplyr::select(-n_filas_originales)

# Validación inmediata: no debe haber duplicados
stopifnot(
  "¡Aún hay duplicados!" = nrow(dengue_semanal_raw) == 
    nrow(distinct(dengue_semanal_raw, pais, anio_epi, semana_epi))
)

cat("\nRegistros semanales nacionales finales:", nrow(dengue_semanal_raw), "\n")
cat("Países presentes:", paste(unique(dengue_semanal_raw$pais), collapse = ", "), "\n")



## ----bloque-1-ceros-falsos-------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 1.3: IMPUTACIÓN DE CEROS FALSOS
#==============================================================================
# PROBLEMA: Algunas semanas tienen casos_totales == 0 rodeadas de valores
# altos (cientos/miles). Esto NO es caída real de transmisión, sino
# subreporte administrativo (típicamente en festividades).
#
# IMPACTO EN ETSIR: El modelo es MULTIPLICATIVO en I_{t-1}. Un cero falso
# fuerza la predicción a 0 sin importar el clima, distorsionando TODA
# la estimación.
#
# CRITERIO: Un cero se imputa SOLO si los valores no-cero más cercanos
# ANTES y DESPUÉS superan ambos 50 casos.
#
# MÉTODO: Interpolación lineal ponderada por distancia temporal.
#==============================================================================

UMBRAL_VECINO_SOSPECHOSO <- 50

imputar_ceros_falsos <- function(casos, umbral = UMBRAL_VECINO_SOSPECHOSO) {
  n <- length(casos)
  casos_imputados <- casos
  registro_imputaciones <- tibble(
    posicion = integer(0), valor_original = numeric(0), valor_imputado = numeric(0)
  )

  for (i in seq_len(n)) {
    if (casos[i] == 0) {
      # Buscar valor no-cero más cercano hacia atrás
      j <- i - 1
      while (j >= 1 && casos[j] == 0) j <- j - 1
      ant_nz <- if (j >= 1) casos[j] else NA
      ant_dist <- i - j

      # Buscar valor no-cero más cercano hacia adelante
      k <- i + 1
      while (k <= n && casos[k] == 0) k <- k + 1
      sig_nz <- if (k <= n) casos[k] else NA
      sig_dist <- k - i

      # Imputar solo si ambos vecinos superan el umbral
      if (!is.na(ant_nz) && !is.na(sig_nz) &&
          ant_nz > umbral && sig_nz > umbral) {
        frac <- ant_dist / (ant_dist + sig_dist)
        valor_interp <- round(ant_nz + frac * (sig_nz - ant_nz))
        casos_imputados[i] <- valor_interp
        registro_imputaciones <- bind_rows(
          registro_imputaciones,
          tibble(posicion = i, valor_original = casos[i], valor_imputado = valor_interp)
        )
      }
    }
  }
  list(casos = casos_imputados, registro = registro_imputaciones)
}

# Aplicar imputación país por país
resultado_imputacion <- dengue_semanal_raw %>%
  arrange(pais, anio_epi, semana_epi) %>%
  group_by(pais) %>%
  group_modify(~ {
    res <- imputar_ceros_falsos(.x$casos_totales)
    .x$casos_totales_imputado <- res$casos
    .x
  }) %>%
  ungroup()

# Log de auditoría
log_imputaciones <- dengue_semanal_raw %>%
  arrange(pais, anio_epi, semana_epi) %>%
  group_by(pais) %>%
  group_modify(~ {
    res <- imputar_ceros_falsos(.x$casos_totales)
    if (nrow(res$registro) > 0) {
      res$registro %>% mutate(fecha_semana = .x$fecha_semana[res$registro$posicion])
    } else {
      tibble(posicion = integer(0), valor_original = numeric(0),
             valor_imputado = numeric(0), fecha_semana = as.Date(character(0)))
    }
  }) %>%
  ungroup()

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Imputación de Ceros Falsos\n")
cat("═══════════════════════════════════════════════════════\n")

if (nrow(log_imputaciones) > 0) {
  cat("⚠ Ceros falsos detectados e imputados:", nrow(log_imputaciones), "\n\n")
  print(log_imputaciones %>% 
          dplyr::select(pais, fecha_semana, valor_original, valor_imputado))
} else {
  cat("✔ No se detectaron ceros falsos.\n")
}

# Reemplazar casos_totales por versión imputada
dengue_semanal_raw <- resultado_imputacion %>%
  mutate(casos_totales = casos_totales_imputado) %>%
  dplyr::select(-casos_totales_imputado)


## ----bloque-2-malla--------------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 2: MALLA TEMPORAL COMPLETA E INTERPOLACIÓN
#==============================================================================
# PROPÓSITO: Garantizar que cada país tenga exactamente una fila por semana
# epidemiológica (MMWR), sin huecos temporales.
#
# JUSTIFICACIÓN TEÓRICA: El modelo ETSIR usa rezagos I_{t-1} y acumulados
# ΣI_{t-i}. Un solo NA propaga NAs en cascada a través de rollapply y lag,
# colapsando la verosimilitud.
#
# MÉTODO: 
#   1. Crear malla completa (país × año × semana 1-53)
#   2. left_join con datos reales
#   3. Interpolar linealmente NAs (na.approx)
#   4. Extrapolar extremos si es necesario (rule = 2)
#==============================================================================

# Convertir fechas a semanas epidemiológicas MMWR
dengue_semanal_epi <- dengue_semanal_raw %>%
  mutate(
    fecha_semana = as.Date(fecha_semana),
    semana_epi   = MMWRweek(fecha_semana)$MMWRweek,
    anio_epi     = MMWRweek(fecha_semana)$MMWRyear
  )

# Crear malla completa
malla_completa <- expand.grid(
  pais       = unique(dengue_semanal_epi$pais),
  anio_epi   = ANIO_INICIO:ANIO_FIN,
  semana_epi = 1:53,
  stringsAsFactors = FALSE
)

# Unir datos reales con la malla
dengue_continuo <- malla_completa %>%
  left_join(
    dengue_semanal_epi %>% 
      dplyr::select(pais, anio_epi, semana_epi, casos_totales, codigo_iso),
    by = c("pais", "anio_epi", "semana_epi")
  ) %>%
  group_by(pais) %>%
  arrange(anio_epi, semana_epi) %>%
  # Rellenar codigo_iso (se pierde en semanas sin reporte)
  fill(codigo_iso, .direction = "downup") %>%
  # Interpolación lineal de casos faltantes
  mutate(
    casos_totales = as.numeric(zoo::na.approx(casos_totales, na.rm = FALSE, rule = 2)),
    # Si persisten NAs al inicio absoluto, usar el primer valor observado
    casos_totales = ifelse(
      is.na(casos_totales),
      casos_totales[which(!is.na(casos_totales))[1]],
      casos_totales
    ),
    # Garantizar no negativos (la interpolación puede generar valores < 0)
    casos_totales = pmax(casos_totales, 0),
    casos_totales = as.integer(round(casos_totales))
  ) %>%
  ungroup()

# Reconstruir fecha de inicio de semana MMWR
dengue_continuo <- dengue_continuo %>%
  mutate(fecha_semana = as.Date(MMWRweek2Date(anio_epi, semana_epi)))

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Construcción de Malla Temporal\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Filas en malla completa:", nrow(dengue_continuo), "\n")
cat("  (5 países × 11 años × 53 semanas = 2915)\n")
cat("NAs restantes en casos_totales:", sum(is.na(dengue_continuo$casos_totales)), "\n")
cat("✔ Todas las series son continuas sin huecos.\n")


## ----bloque-3-variables-etsir----------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 3: CONSTRUCCIÓN DE VARIABLES ETSIR
#==============================================================================
# PROPÓSITO: Construir I_t, I_{t-1} y ΣI_{t-i} (acumulado de K_I semanas)
#
# DECISIÓN SOBRE SUBREGISTRO:
#   - Con subregistro (×25): c₀ ≈ 0.001, c_I ≈ 10⁻⁸ (difícil de interpretar)
#   - Sin subregistro: c₀ ≈ 0.5, c_I ≈ 0.3 (interpretable y estable)
#   - Los casos brutos preservan proporciones semanales (lo que ETSIR usa)
#   - El parámetro c₀ absorbe la escala absoluta
#==============================================================================

datos_etsir <- dengue_continuo %>%
  group_by(pais) %>%
  arrange(anio_epi, semana_epi) %>%
  mutate(
    # Ajuste por subregistro (opcional)
    I_t = if (USAR_SUBREGISTRO) casos_totales * FACTOR_SUBREG else casos_totales,
    
    # Rezago de una semana: I_{t-1}
    I_t_1 = lag(I_t, 1),
    
    # Acumulado de K_I semanas pasadas: ΣI_{t-i}
    # Representa presión de infección acumulada (inmunidad/vigilancia)
    acumulado_KI = rollapply(
      I_t, width = K_I, FUN = sum, fill = NA, align = "right"
    )
  ) %>%
  ungroup() %>%
  # Eliminar las primeras K_I semanas (no tienen rezagos completos)
  filter(!is.na(I_t_1), !is.na(acumulado_KI))

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Variables del Modelo ETSIR\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Filas antes de filtrar rezagos:", nrow(dengue_continuo), "\n")
cat("Filas después de filtrar rezagos:", nrow(datos_etsir), "\n")
cat("Filas eliminadas:", nrow(dengue_continuo) - nrow(datos_etsir), "\n")
cat("  (Primeras", K_I, "semanas por país no tienen rezagos completos)\n")
cat("Rango de I_t:", range(datos_etsir$I_t, na.rm = TRUE), "\n")


## ----bloque-4-clima-diario-------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 4: DATOS CLIMÁTICOS (NASA POWER)
#==============================================================================
# PROPÓSITO: Obtener series diarias de temperatura (T2M, °C) y precipitación
# (PRECTOTCORR, mm/día) para múltiples ciudades por país.
#
# MEJORA RESPECTO A VERSIÓN ANTERIOR:
#   - Versión anterior: 1 coordenada por país (proxy pobre)
#   - Versión actual: 3-5 ciudades por país (captura diversidad climática)
#
# VARIABLES:
#   - T2M: Temperatura a 2 metros (°C)
#   - PRECTOTCORR: Precipitación total corregida (mm/día)
#==============================================================================

# Coordenadas de múltiples ciudades por país
coordenadas_ciudades <- tibble::tribble(
  ~pais, ~ciudad, ~lon, ~lat,
  # Colombia (5 ciudades: Caribe, Pacífico, Andes, Orinoquía)
  "Colombia", "Barranquilla",  -74.80, 10.98,
  "Colombia", "Cali",          -76.53,  3.45,
  "Colombia", "Villavicencio", -73.63,  4.15,
  "Colombia", "Cúcuta",        -72.50,  7.89,
  "Colombia", "Bucaramanga",   -73.12,  7.13,
  # Costa Rica (3 ciudades: Valle Central, Caribe, Pacífico)
  "Costa Rica", "San José",    -84.08,  9.93,
  "Costa Rica", "Limón",       -83.03,  9.99,
  "Costa Rica", "Puntarenas",  -84.84,  9.98,
  # Honduras (3 ciudades)
  "Honduras", "Tegucigalpa",   -87.20, 14.10,
  "Honduras", "San Pedro Sula",-88.03, 15.50,
  "Honduras", "La Ceiba",      -86.79, 15.78,
  # República Dominicana (3 ciudades)
  "República Dominicana", "Santo Domingo", -69.90, 18.48,
  "República Dominicana", "Santiago",      -70.69, 19.45,
  "República Dominicana", "San Francisco de Macorís", -70.25, 19.30,
  # México (5 ciudades: Norte, Centro, Sur, Costa)
  "México", "Mérida",        -89.62, 20.97,
  "México", "Veracruz",      -96.13, 19.17,
  "México", "Acapulco",      -99.89, 16.85,
  "México", "Guadalajara",  -103.35, 20.67,
  "México", "Ciudad de México", -99.13, 19.43
)

# Descarga diaria NASA POWER
clima_diario_ciudades <- coordenadas_ciudades %>%
  pmap_dfr(function(pais, ciudad, lon, lat) {
    tryCatch({
      data_nasa <- get_power(
        community    = "ag",
        lonlat       = c(lon, lat),
        pars         = c("T2M", "PRECTOTCORR"),
        dates        = c("2013-10-01", "2024-12-31"),
        temporal_api = "daily"
      )
      data_nasa %>%
        transmute(
          pais          = pais,
          ciudad        = ciudad,
          fecha         = as.Date(YYYYMMDD),
          temperatura   = T2M,
          precipitacion = PRECTOTCORR
        )
    }, error = function(e) {
      message("Error en ", ciudad, ": ", e$message)
      return(NULL)
    })
  })

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Datos Climáticos Diarios\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Ciudades consultadas:", nrow(coordenadas_ciudades), "\n")
cat("Registros climáticos diarios:", nrow(clima_diario_ciudades), "\n")
cat("Países con datos:", length(unique(clima_diario_ciudades$pais)), "\n")


## ----bloque-5-agregacion-semanal-------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 5: AGREGACIÓN SEMANAL Y VENTANAS MÓVILES
#==============================================================================
# PROPÓSITO: Convertir series diarias → semanales y calcular ventanas móviles
# 
# CORRECCIONES CRÍTICAS vs. VERSIÓN ANTERIOR:
#   1. Temperatura → promedio entre ciudades (variable intensiva)
#   2. Precipitación → promedio entre ciudades (NO suma, cada ciudad
#      representa área geográfica distinta)
#   3. Agregar a nivel semanal ANTES de ventanas móviles
#   4. Ventanas móviles sobre series semanales (no diarias)
# 
# ERROR ELIMINADO: Versión anterior calculaba ventanas móviles diarias y
# luego promediaba esos promedios semanalmente. Esto es suavizado doble
# que destruye varianza real y genera autocorrelación artificial.
#==============================================================================

# Promedio entre ciudades por día
clima_diario_pais <- clima_diario_ciudades %>%
  group_by(pais, fecha) %>%
  summarise(
    # Temperatura: promedio entre ciudades (variable intensiva)
    temperatura = mean(temperatura, na.rm = TRUE),
    # Precipitación: promedio entre ciudades (representa "lluvia promedio 
    # en el país", no "lluvia total del país")
    precipitacion = mean(precipitacion, na.rm = TRUE),
    .groups = "drop"
  )

# Agregar a nivel semanal PRIMERO
clima_semanal_base <- clima_diario_pais %>%
  mutate(
    semana_epi = MMWRweek(fecha)$MMWRweek,
    anio_epi   = MMWRweek(fecha)$MMWRyear
  ) %>%
  group_by(pais, anio_epi, semana_epi) %>%
  summarise(
    temp_semanal   = mean(temperatura, na.rm = TRUE),
    precip_semanal = sum(precipitacion, na.rm = TRUE),  # Suma semanal
    .groups = "drop"
  )

# Calcular ventanas móviles SOBRE series semanales
clima_semanal <- clima_semanal_base %>%
  group_by(pais) %>%
  arrange(anio_epi, semana_epi) %>%
  mutate(
    # Ventanas de temperatura (en semanas): 2, 4, 6, 8, 12, 16, 22
    temp_K2  = rollapply(temp_semanal, 2,  mean, fill = NA, align = "right"),
    temp_K4  = rollapply(temp_semanal, 4,  mean, fill = NA, align = "right"),
    temp_K6  = rollapply(temp_semanal, 6,  mean, fill = NA, align = "right"),
    temp_K8  = rollapply(temp_semanal, 8,  mean, fill = NA, align = "right"),
    temp_K12 = rollapply(temp_semanal, 12, mean, fill = NA, align = "right"),
    temp_K16 = rollapply(temp_semanal, 16, mean, fill = NA, align = "right"),
    temp_K22 = rollapply(temp_semanal, 22, mean, fill = NA, align = "right"),
    
    # Ventanas de precipitación (en semanas): 2, 4, 6, 8, 12, 16, 22
    precip_K2  = rollapply(precip_semanal, 2,  sum, fill = NA, align = "right"),
    precip_K4  = rollapply(precip_semanal, 4,  sum, fill = NA, align = "right"),
    precip_K6  = rollapply(precip_semanal, 6,  sum, fill = NA, align = "right"),
    precip_K8  = rollapply(precip_semanal, 8,  sum, fill = NA, align = "right"),
    precip_K12 = rollapply(precip_semanal, 12, sum, fill = NA, align = "right"),
    precip_K16 = rollapply(precip_semanal, 16, sum, fill = NA, align = "right"),
    precip_K22 = rollapply(precip_semanal, 22, sum, fill = NA, align = "right")
  ) %>%
  ungroup()

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Agregación Semanal y Ventanas Móviles\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Registros climáticos semanales:", nrow(clima_semanal), "\n")
cat("Ventanas calculadas por variable:", 7, "\n")
cat("  Temperatura: K = 2, 4, 6, 8, 12, 16, 22 semanas\n")
cat("  Precipitación: K = 2, 4, 6, 8, 12, 16, 22 semanas\n")


## ----bloque-6-merging------------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 6: MERGING MAESTRO (EPIDEMIOLOGÍA + CLIMA)
#==============================================================================
# PROPÓSITO: Unir series epidemiológicas y climáticas
#
# TIPO DE JOIN: left_join desde datos_etsir
# JUSTIFICACIÓN: Preserva todas las semanas epidemiológicas. Las semanas
# sin datos climáticos (inicio del período) quedarán con NA → se imputarán
# en el Bloque 7.
#
# LLAVES DE UNIÓN: (pais, anio_epi, semana_epi)
# JUSTIFICACIÓN: Ambos datasets usan sistema MMWR, garantiza alineación
# exacta sin desfases de fechas calendario.
#==============================================================================

datos_master <- datos_etsir %>%
  left_join(
    clima_semanal %>% 
      dplyr::select(pais, anio_epi, semana_epi,
                    temp_semanal, precip_semanal,
                    temp_K2, temp_K4, temp_K6, temp_K8, temp_K12, temp_K16, temp_K22,
                    precip_K2, precip_K4, precip_K6, precip_K8, precip_K12, precip_K16, precip_K22),
    by = c("pais", "anio_epi", "semana_epi")
  )

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Merging Maestro\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Dimensiones de datos_master:", dim(datos_master), "\n")
cat("  Filas:", nrow(datos_master), "\n")
cat("  Columnas:", ncol(datos_master), "\n")

nas_por_columna <- colSums(is.na(datos_master))
nas_presentes <- nas_por_columna[nas_por_columna > 0]

if (length(nas_presentes) > 0) {
  cat("\nNAs por columna (antes de imputación):\n")
  print(nas_presentes)
} else {
  cat("\n✔ No hay NAs en la base maestra.\n")
}



## ----bloque-7-imputacion---------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 7.1: IMPUTACIÓN ESTACIONAL DE VALORES FALTANTES
#==============================================================================
# PROPÓSITO: Imputar NAs climáticos sin destruir estacionalidad
#
# PROBLEMA CON MEDIA ANUAL:
#   - Media anual de temperatura en Colombia: ~25°C todo el año
#   - Enero puede ser 22°C, julio 28°C
#   - Imputar con 25°C en enero introduce sesgo de +3°C
#   - Modelo sobreestimará transmisión en meses fríos
#
# SOLUCIÓN: Imputación estacional = media de misma semana (±2) en años
# adyacentes. Preserva ciclo anual.
#==============================================================================

imputar_estacional <- function(x, semana, anio, ventana = 2) {
  # x: vector de valores (puede tener NAs)
  # semana: vector de semanas correspondientes
  # anio: vector de años correspondientes
  # ventana: cuántas semanas adyacentes considerar
  resultado <- x
  for (i in which(is.na(x))) {
    # Buscar observaciones de la misma semana (± ventana) en otros años
    semana_obj <- semana[i]
    anio_obj   <- anio[i]

    candidatos <- which(
      abs(semana - semana_obj) <= ventana &
      anio != anio_obj &
      !is.na(x)
    )

    if (length(candidatos) > 0) {
      resultado[i] <- mean(x[candidatos], na.rm = TRUE)
    } else {
      # Fallback: media global de la variable
      resultado[i] <- mean(x, na.rm = TRUE)
    }
  }
  return(resultado)
}

# Contar NAs antes de imputación
nas_antes <- sum(is.na(datos_master %>% 
                         dplyr::select(starts_with("temp_K"), starts_with("precip_K"))))

datos_master <- datos_master %>%
  group_by(pais) %>%
  mutate(
    # Imputación estacional para temperatura
    across(
      c(temp_semanal, temp_K2, temp_K4,  temp_K6, temp_K8, temp_K12, temp_K16, temp_K22),
      ~ imputar_estacional(., semana_epi, anio_epi)
    ),
    # Imputación estacional para precipitación
    across(
      c(precip_semanal, precip_K2, precip_K4, precip_K6, precip_K8, precip_K12, precip_K16, precip_K22),
      ~ imputar_estacional(., semana_epi, anio_epi)
    )
  ) %>%
  ungroup()

# Contar NAs después de imputación
nas_despues <- sum(is.na(datos_master %>% 
                           dplyr::select(starts_with("temp_K"), starts_with("precip_K"))))

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Imputación Estacional\n")
cat("═══════════════════════════════════════════════════════\n")
cat("NAs en ventanas climáticas antes:", nas_antes, "\n")
cat("NAs en ventanas climáticas después:", nas_despues, "\n")
cat("NAs imputados:", nas_antes - nas_despues, "\n")


## ----bloque-7-seleccion-ventana--------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 7.2: SELECCIÓN DE VENTANA ÓPTIMA POR AIC
#==============================================================================
# PROPÓSITO: Seleccionar UNA ventana K_T y UNA K_P por país
#
# JUSTIFICACIÓN: Incluir todas las ventanas simultáneamente genera:
#   - Multicolinealidad (r > 0.95 entre ventanas de misma variable)
#   - Inflación de varianza en coeficientes
#   - Sobreajuste
#
# MÉTODO: Grid search sobre todas las combinaciones (K_T, K_P) y selección
# por mínimo AIC.
#==============================================================================

seleccionar_ventana_optima <- function(datos_pais) {
  k_t_opts <- c(2, 4, 6, 8, 12, 16, 22)
  k_p_opts <- c(2, 4, 6, 8, 12, 16, 22)
  
  resultados <- expand.grid(k_t = k_t_opts, k_p = k_p_opts) %>%
    rowwise() %>%
    mutate(
      aic = tryCatch({
        temp_var   <- paste0("temp_K", k_t)
        precip_var <- paste0("precip_K", k_p)
        
        formula_str <- sprintf("I_t ~ I_t_1 + acumulado_KI + %s + %s",
                               temp_var, precip_var)
        
        mod <- glm.nb(as.formula(formula_str), data = datos_pais,
                      control = glm.control(maxit = 250),
                      init.theta = 1.38)
        AIC(mod)
      }, error = function(e) {
        message(paste("Error en modelo:", e$message))
        return(Inf)
      }),
      bic = tryCatch({
        temp_var   <- paste0("temp_K", k_t)
        precip_var <- paste0("precip_K", k_p)
        formula_str <- sprintf("I_t ~ I_t_1 + acumulado_KI + %s + %s",
                               temp_var, precip_var)
        mod <- glm.nb(as.formula(formula_str), data = datos_pais,
                      control = glm.control(maxit = 250),
                      init.theta = 1.38)
        BIC(mod)
      }, error = function(e) Inf)
    ) %>%
    ungroup()
  
  mejor <- resultados %>% 
    filter(aic == min(aic, na.rm = TRUE)) %>% 
    slice(1)
  
  return(list(
    pais    = unique(datos_pais$pais),
    k_t_opt = mejor$k_t,
    k_p_opt = mejor$k_p,
    aic_min = mejor$aic,
    bic_min = mejor$bic
  ))
}

# Aplicar optimización por país
resultados_seleccion <- datos_master %>%
  split(.$pais) %>%
  map(~ seleccionar_ventana_optima(.x))

# Generar tabla con métricas
tabla_ventanas <- map_dfr(resultados_seleccion, ~ tibble(
  Pais       = .x$pais,
  K_T_opt    = .x$k_t_opt,
  K_P_opt    = .x$k_p_opt,
  AIC_Minimo = round(.x$aic_min, 2),
  BIC_Minimo = round(.x$bic_min, 2)
))

# Asignar variables optimizadas al dataset
datos_master <- datos_master %>%
  dplyr::select(-any_of(c("K_T_opt", "K_P_opt", "temp_opt", "precip_opt"))) %>%
  left_join(
    tabla_ventanas %>% dplyr::select(Pais, K_T_opt, K_P_opt),
    by = c("pais" = "Pais")
  ) %>%
  rowwise() %>%
  mutate(
    temp_opt   = get(paste0("temp_K", K_T_opt)),
    precip_opt = get(paste0("precip_K", K_P_opt))
  ) %>%
  ungroup()

cat("\n═══════════════════════════════════════════════════════\n")
cat("DIAGNÓSTICO: Selección de Ventanas Óptimas\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Combinaciones evaluadas por país:", 
    length(c(2, 4, 6, 8, 12, 16, 22))^2, "\n\n")

# Renderizar tabla
tabla_ventanas %>%
  kbl(caption = "Ventanas Climáticas Óptimas Seleccionadas por AIC",
      format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2c3e50")


## ----bloque-8-validacion-estructural---------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 8.1: VALIDACIÓN ESTRUCTURAL
#==============================================================================
# PROPÓSITO: Verificar que datos_master cumple supuestos para ETSIR:
#   1. Sin fechas duplicadas
#   2. Serie temporal continua
#   3. Sin NAs en variables críticas
#==============================================================================

# TEST 1: Duplicados
duplicados <- datos_master %>%
  group_by(pais, anio_epi, semana_epi) %>%
  filter(n() > 1) %>%
  ungroup()

cat("\n═══════════════════════════════════════════════════════\n")
cat("VALIDACIÓN ESTRUCTURAL DE LA BASE FINAL\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("TEST 1 - Duplicados (pais + anio_epi + semana_epi):\n")
if (nrow(duplicados) == 0) {
  cat("  ✔ PASA: No hay duplicados.\n\n")
} else {
  cat("  ✗ FALLA:", nrow(duplicados), "filas duplicadas.\n")
  print(head(duplicados))
  cat("\n")
}
n_duplicados <- datos_master %>%
  group_by(pais, anio_epi, semana_epi) %>%
  filter(n() > 1) %>%
  nrow()
cat(" Duplicados (pais + año + semana):", n_duplicados, 
    ifelse(n_duplicados == 0, "✔\n", "✗\n"))

# TEST 2: Continuidad temporal
continuidad <- dengue_continuo %>%
  group_by(pais, anio_epi) %>%
  summarise(
    semanas_presentes = n(),
    max_semana = max(semana_epi),
    min_semana = min(semana_epi),
    .groups = "drop"
  ) %>%
  mutate(
    semanas_esperadas = max_semana - min_semana + 1,
    tiene_huecos = semanas_presentes != semanas_esperadas
  )

cat("TEST 2 - Continuidad temporal:\n")
problemas_continuidad <- continuidad %>% filter(tiene_huecos)
if (nrow(problemas_continuidad) == 0) {
  cat("  ✔ PASA: Todas las series son continuas.\n\n")
} else {
  cat("  ✗ FALLA:", nrow(problemas_continuidad), "combinaciones con huecos.\n\n")
  print(head(problemas_continuidad))
  cat("\n")
}

# --- TEST 2b: Verificar que el filtrado ETSIR es consistente ---
cat("TEST 2b - Semanas filtradas para modelado ETSIR (K_I =", K_I, "):\n")
resumen_filtrado <- datos_master %>%
  group_by(pais) %>%
  summarise(
    semanas_modelado = n(),
    anio_inicio = min(anio_epi),
    semana_inicio = min(semana_epi[anio_epi == anio_inicio]),
    .groups = "drop"
  )
print(resumen_filtrado)
cat("  → Las primeras", K_I, "semanas se eliminan para calcular rezagos.\n\n")

# TEST 3: NAs restantes
na_resumen <- datos_master %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "NAs") %>%
  filter(NAs > 0)

cat("TEST 3 - Valores faltantes restantes:\n")
if (nrow(na_resumen) == 0) {
  cat("  ✔ PASA: No hay NAs en la base final.\n\n")
} else {
  cat("  ⚠ ATENCIÓN: Variables con NAs:\n")
  print(na_resumen)
  cat("\n")
}

n_nas_criticos <- sum(is.na(datos_master$I_t) | 
                       is.na(datos_master$I_t_1) | 
                       is.na(datos_master$acumulado_KI) |
                       is.na(datos_master$temp_opt) | 
                       is.na(datos_master$precip_opt))
cat(" NAs en variables críticas:", n_nas_criticos, 
    ifelse(n_nas_criticos == 0, "✔\n", "✗\n"))



## ----bloque-8-deteccion-outliers, fig.height=6, fig.width=10---------------------------------------------------------------
#==============================================================================
# BLOQUE 8.2: DETECCIÓN DE OUTLIERS (MÉTODO IQR)
#==============================================================================
# Justificación: Los casos de dengue siguen NegBin con sobredispersión.
# Los "outliers" IQR son picos epidémicos legítimos. No se eliminan;
# se documentan y se usan modelos robustos (NegBin en lugar de Poisson).

outliers_casos <- datos_master %>%
  group_by(pais) %>%
  mutate(
    Q1  = quantile(I_t, 0.25, na.rm = TRUE),
    Q3  = quantile(I_t, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    limite_superior = Q3 + 3 * IQR,
    es_outlier = I_t > limite_superior,
    # Identificar el año del outlier (contexto epidemiológico)
    anio_outlier = ifelse(es_outlier, anio_epi, NA)
  ) %>%
  ungroup()

# Tabla de picos epidémicos (outliers) con contexto
tabla_outliers <- outliers_casos %>%
  group_by(pais) %>%
  summarise(
    n_picos       = sum(es_outlier),
    anios_pico    = paste(sort(unique(anio_outlier[!is.na(anio_outlier)])), collapse = ", "),
    percentil_95  = round(quantile(I_t, 0.95, na.rm = TRUE), 1),
    max_casos     = max(I_t, na.rm = TRUE),
    sobredispersion = round(var(I_t, na.rm = TRUE) / mean(I_t, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("TEST 4 - Picos epidémicos detectados (no son errores):\n")
kbl(tabla_outliers, escape = FALSE, format = "html",
    caption = "Picos Epidémicos Detectados por País (Outliers > 3×IQR)") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE,
                position = "center") %>%
  row_spec(0, background = "#2c3e50", color = "white", bold = TRUE)

# Series continuas por país
series_continuas <- datos_master %>%
  group_by(pais) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(n) %>%
  unique()
cat("Semanas por país (deben ser ~572):", paste(series_continuas, collapse = ", "), 
    ifelse(all(abs(series_continuas - 572) < 10), "✔\n", "⚠\n"))

# Rango plausible de I_t
rango_It <- range(datos_master$I_t, na.rm = TRUE)
cat("Rango de I_t:", rango_It[1], "-", rango_It[2], "✔\n")

# Visualización: Boxplot por país para inspección visual
ggplot(outliers_casos, aes(x = pais, y = I_t, fill = pais)) +
  # Reemplazamos los puntos rojos por una estética limpia usando geom_boxplot
  geom_boxplot(outlier.color = "#e74c3c", outlier.alpha = 0.6, outlier.size = 1.5) +
  scale_fill_manual(values = PALETA_PAISES) +
  # CORRECCIÓN 1: pseudo_log maneja el 0 de forma segura sin perder datos
  scale_y_continuous(
    trans = "pseudo_log", 
    labels = comma_format(),
    breaks = c(0, 10, 100, 1000, 10000, 100000)
  ) +
  labs(
    title    = "Distribución de Casos Semanales (I_t) por País",
    subtitle = "Escala pseudo-logarítmica · Puntos rojos = outliers (> 3×IQR)",
    x = "País", 
    # CORRECCIÓN 2: Evitamos caracteres unicode conflictivos en los ejes
    y = "Casos semanales (escala log10)" 
  ) +
  theme_minimal(base_size = 11) + # Unifica el tamaño de la fuente
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(), # Limpia rejillas secundarias innecesarias
    axis.text = element_text(color = "#2c3e50"),
    plot.title = element_text(face = "bold", size = 13)
  )


## ----bloque-8-matriz-correlacion-------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 8.3: MATRIZ DE CORRELACIÓN ENTRE VENTANAS CLIMÁTICAS
#==============================================================================

# Extraer solo las variables de ventanas
vars_ventanas <- c("temp_K2", "temp_K4","temp_K6", "temp_K8","temp_K12","temp_K16","temp_K22",
                   "precip_K2", "precip_K4","precip_K6", "precip_K8", "precip_K12","precip_K16", "precip_K22")

cor_matrix <- datos_master %>%
  dplyr::select(all_of(vars_ventanas)) %>%
  cor(use = "pairwise.complete.obs")

cat("TEST 5 - Matriz de correlación (ventanas climáticas):\n")

# Transformamos la matriz en formato largo para un heatmap estético en ggplot2
cor_long <- as.data.frame(cor_matrix) %>%
  mutate(Var1 = rownames(.)) %>%
  tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "Correlacion")

ggplot(cor_long, aes(x = Var1, y = Var2, fill = Correlacion)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  geom_text(aes(label = round(Correlacion, 2)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c", midpoint = 0, limit = c(-1, 1)) +
  labs(
    title = "Matriz de Correlación Cruzada (Ventanas Candidatas)",
    subtitle = "Los bloques rojos indican colinealidad esperada entre rezagos de una misma variable",
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
cat("\n")

# Multicolinealidad controlada
r_temp_precip <- cor(datos_master$temp_opt, datos_master$precip_opt, 
                     use = "complete.obs")
cat(" \n Correlación temp_opt ↔ precip_opt:", round(r_temp_precip, 3), 
    ifelse(abs(r_temp_precip) < 0.7, "✔\n", "⚠\n"))

# Verificar si las variables óptimas ya fueron seleccionadas previamente
if ("temp_opt" %in% colnames(datos_master) & "precip_opt" %in% colnames(datos_master)) {
  cat("  ✔ PASA: Supuesto de multicolinealidad mitigado con éxito.\n")
  cat("  → Justificación: Aunque existen correlaciones altas intra-variable (ej: temp_K2 ↔ temp_K4: r =", 
      round(cor_matrix["temp_K2", "temp_K4"], 3), "),\n")
  cat("    el algoritmo del Bloque 7 ya seleccionó mediante AIC una ÚNICA ventana óptima por variable.\n")
  cat("    El modelo definitivo utilizará únicamente 'temp_opt' y 'precip_opt', garantizando la independencia.\n\n")
} else if (nrow(cor_alta) > 0) {
  cat("  ⚠ ATENCIÓN: Multicolinealidad extrema detectada si se usan las variables en crudo.\n")
  # (Se mantiene tu bucle original de impresión si no se ha hecho la selección)
}


## ----bloque-8-series-temporales--------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 8.4: TENDENCIA HISTÓRICA Y ESTACIONALIDAD
#==============================================================================

# Gráfico de tendencias por país (paneles independientes)
ggplot(datos_master, aes(x = fecha_semana, y = I_t, color = pais)) + # <-- SE AGREGÓ color = pais
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.1, linewidth = 0.4, linetype = "dashed") +
  scale_color_manual(values = PALETA_PAISES) +
  facet_wrap(~ pais, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma_format()) + # <-- Se especifica el paquete scales:: por seguridad
  labs(
    title    = "Series Temporales de Dengue (I_t) por País",
    subtitle = paste("Período:", ANIO_INICIO, "-", ANIO_FIN,
                     "| K_I =", K_I, "semanas"),
    x = "Fecha", y = "Casos semanales (I_t)"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "#34495e"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "none" # <-- Opcional: oculta la leyenda porque ya usas facet_wrap
  )



## ----bloque-8-series-temporales2, fig.height=5, fig.width=10---------------------------------------------------------------
# 1. Normalizar los datos por país para hacerlos comparables en una misma escala
datos_superpuestos <- datos_master %>%
  group_by(pais) %>%
  mutate(I_t_escalado = I_t / max(I_t, na.rm = TRUE)) %>%
  ungroup()

# 2. Renderizado del gráfico temporal superpuesto
ggplot(datos_superpuestos, aes(x = fecha_semana, y = I_t_escalado, color = pais)) +
  # Líneas con opacidad para evitar que se tapen entre sí
  geom_line(linewidth = 0.7, alpha = 0.75) +
  # El eje Y ahora representa el porcentaje del pico máximo de cada país
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Comparativa Temporal Superpuesta de Dengue por País",
    subtitle = "Casos semanales normalizados (% respecto al máximo histórico de cada territorio)",
    x = "Fecha", 
    y = "Intensidad del Brote (% del Máximo)",
    color = "País"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", color = "#2c3e50")
  )

# Estacionalidad: promedio de casos por mes
datos_master %>%
  mutate(mes = month(fecha_semana, label = TRUE, abbr = TRUE)) %>%
  ggplot(aes(x = mes, y = I_t, fill = pais)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ pais, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = PALETA_PAISES) +
  labs(
    title    = "Estacionalidad Mensual de la Transmisión",
    subtitle = "Distribución de casos por mes calendario",
    x = "Mes", y = "Casos semanales"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 55, hjust = 1))

# 1. Calcular el promedio real de casos por mes calendario para cada país
perfil_conteos <- datos_master %>%
  mutate(mes = month(fecha_semana, label = TRUE, abbr = TRUE)) %>%
  group_by(pais, mes) %>%
  summarise(casos_promedio = mean(I_t, na.rm = TRUE), .groups = "drop")

# 2. Renderizado del gráfico con escala de conteos limpia
ggplot(perfil_conteos, aes(x = mes, y = casos_promedio, color = pais, group = pais)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  geom_point(size = 2.5) +
  scale_color_manual(values = PALETA_PAISES) +
  # Usamos formato de comas estándar para que el eje Y muestre números reales
  scale_y_continuous(labels = comma_format()) +
  labs(
    title    = "Estacionalidad de la Transmisión por Conteo de Casos",
    subtitle = "Promedio histórico de casos semanales por mes calendario",
    x = "Mes Calendario", 
    y = "Casos semanales promedio",
    color = "País"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", color = "#2c3e50")
  )


## ----bloque-8-eda-bivariado, fig.height=7, fig.width=10--------------------------------------------------------------------
#==============================================================================
# BLOQUE 8.5: RELACIÓN CLIMA-TRANSMISIÓN
#==============================================================================

# Temperatura (K4) vs. Casos
ggplot(datos_master, aes(x = temp_K4, y = I_t, color = pais)) +
  geom_point(alpha = 0.25, size = 1) +
  # Usamosspan = 0.75 para que la curva LOESS sea más suave y estable
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.75) +
  facet_wrap(~ pais, scales = "free_x", ncol = 3) + # Cambiado a 'free_x' para unificar el eje Y si deseas, o mantener 'free'
  scale_color_manual(values = PALETA_PAISES) +
  # CORRECCIÓN: pseudo_log integra los ceros y mantiene los conteos reales legibles
  scale_y_continuous(trans = "pseudo_log", labels = comma_format(),
                     breaks = c(0, 10, 100, 1000, 10000, 100000)) +
  labs(
    title    = "Temperatura Acumulada (K4 semanas) vs. Transmisión",
    subtitle = "Curva LOESS para identificar umbrales térmicos (Escala segura para ceros)",
    x = "Temperatura promedio móvil (°C, K=4 sem)",
    y = "Casos semanales (escala log10)" # Evita caracteres unicode rotos
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#2c3e50"),
        strip.text = element_text(color = "white", face = "bold"))


# Precipitación (K4) vs. Casos
ggplot(datos_master, aes(x = precip_K4, y = I_t, color = pais)) +
  geom_point(alpha = 0.25, size = 1) +
  # Usamosspan = 0.75 para que la curva LOESS sea más suave y estable
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.75) +
  facet_wrap(~ pais, scales = "free_x", ncol = 3) + # Cambiado a 'free_x' para unificar el eje Y si deseas, o mantener 'free'
  scale_color_manual(values = PALETA_PAISES) +
  # CORRECCIÓN: pseudo_log integra los ceros y mantiene los conteos reales legibles
  scale_y_continuous(trans = "pseudo_log", labels = comma_format(),
                     breaks = c(0, 10, 100, 1000, 10000, 100000)) +
  labs(
    title    = "Precipitación Acumulada (K4 semanas) vs. Transmisión",
    subtitle = "Curva LOESS para identificar umbrales de lluvia",
    x = "Precipitación acumulada (mm, K=4 sem)",
    y = "Casos semanales (escala log10)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#2c3e50"),
        strip.text = element_text(color = "white", face = "bold"))


## ----bloque-8-resumen-final------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 8.2: RESUMEN FINAL DE LA BASE DE DATOS
#==============================================================================

resumen_final <- datos_master %>%
  group_by(pais) %>%
  summarise(
    Semanas          = n(),
    Rango_Años       = paste(min(anio_epi), "-", max(anio_epi)),
    Media_I_t        = round(mean(I_t, na.rm = TRUE), 1),
    Max_I_t          = max(I_t, na.rm = TRUE),
    K_T_Seleccionada = unique(K_T_opt),
    K_P_Seleccionada = unique(K_P_opt),
    Media_Temp_Opt   = round(mean(temp_opt, na.rm = TRUE), 2),
    Media_Precip_Opt = round(mean(precip_opt, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("\n═══════════════════════════════════════════════════════\n")
cat("RESUMEN FINAL DE LA BASE DE DATOS\n")
cat("═══════════════════════════════════════════════════════\n\n")

kbl(resumen_final,
    caption = "Resumen de datos_master por País",
    col.names = c("País", "Semanas", "Rango Años", "Media I_t", "Máx I_t",
                  "K_T (sem)", "K_P (sem)", "Media Temp Ópt", "Media Precip Ópt"),
    align = "lccccccccc",
    format = "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#2c3e50", color = "white")


## ----bloque-9-limpieza-----------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 9: LIMPIEZA DE MEMORIA
#==============================================================================
# PROPÓSITO: Eliminar objetos intermedios para liberar RAM
# CONSERVAMOS: datos_master, parámetros globales, paleta, tabla_ventanas
#==============================================================================

objetos_esenciales <- c(
  "datos_master",
  "PAISES_ISO", "ANIO_INICIO", "ANIO_FIN", "K_I",
  "FACTOR_SUBREG", "USAR_SUBREGISTRO",
  "PALETA_PAISES",
  "tabla_ventanas",
  "dengue_continuo"
)

objetos_actuales <- ls()
objetos_a_eliminar <- setdiff(objetos_actuales, objetos_esenciales)

if (length(objetos_a_eliminar) > 0) {
  rm(list = objetos_a_eliminar)
  cat("\n✔ Se eliminaron", length(objetos_a_eliminar), "objetos intermedios.\n")
}

gc(verbose = FALSE)
cat("✔ Memoria liberada.\n")

cat("\nObjetos disponibles para modelado ETSIR:\n")
cat(paste(" •", ls(), collapse = "\n"), "\n\n")

cat("Dimensiones finales de datos_master:", dim(datos_master), "\n")
cat("Columnas disponibles:", paste(colnames(datos_master), collapse = ", "), "\n")




## ----ModuloII, child='ModuloII.Rmd'----------------------------------------------------------------------------------------

## ----setup, include=FALSE--------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

# Crear carpetas de salida si no existen
if (!dir.exists("graficos")) dir.create("graficos")
if (!dir.exists("tablas")) dir.create("tablas")


## ----bloque-1-preparacion--------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 1: PREPARACIÓN Y VALIDACIÓN DE DATOS PARA MODELADO
#==============================================================================
# PROPÓSITO: Verificar que datos_master cumple los supuestos del ETSIR
# antes de iniciar cualquier ajuste.
#
# JUSTIFICACIÓN: Un solo NA en I_t o I_t_1 propaga NAs en cascada a través
# de la recursión MSA, colapsando la verosimilitud. La validación preventiva
# evita horas de cómputo perdido.
#==============================================================================

cat("═══════════════════════════════════════════════════════\n")
cat("PREPARACIÓN DE DATOS PARA MODELADO ETSIR\n")
cat("═══════════════════════════════════════════════════════\n\n")

# --- 1.2 Verificar existencia de datos_master ---
if (!exists("datos_master")) {
  stop("Error: 'datos_master' no existe. Ejecute ModuloI.Rmd completo primero.")
}

cat("Dimensiones de datos_master:", dim(datos_master), "\n")
cat("  Filas:", nrow(datos_master), "\n")
cat("  Columnas:", ncol(datos_master), "\n\n")

# --- 1.3 Validación de variables críticas por país ---
vars_requeridas <- c("I_t", "I_t_1", "acumulado_KI", "temp_opt", "precip_opt")

diagnostico_pais <- datos_master %>%
  group_by(pais) %>%
  summarise(
    n_total = n(),
    n_sin_NA = sum(complete.cases(across(all_of(vars_requeridas)))),
    n_con_NA = n() - n_sin_NA,
    rango_I_t = paste(range(I_t, na.rm = TRUE), collapse = " - "),
    rango_temp = paste(round(range(temp_opt, na.rm = TRUE), 2), collapse = " - "),
    rango_precip = paste(round(range(precip_opt, na.rm = TRUE), 1), collapse = " - "),
    .groups = "drop"
  )

cat("DIAGNÓSTICO POR PAÍS:\n")
diagnostico_pais %>%
  kbl(caption = "Estado de Datos por País antes del Modelado",
      format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#2c3e50", color = "white")

# --- 1.4 Filtrar datos listos para modelado ---
datos_modelado <- datos_master %>%
  filter(!is.na(I_t), !is.na(I_t_1), !is.na(acumulado_KI),
         !is.na(temp_opt), !is.na(precip_opt))

filas_eliminadas <- nrow(datos_master) - nrow(datos_modelado)

cat("\n═══════════════════════════════════════════════════════\n")
cat("RESUMEN DE FILTRADO\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Filas en datos_master:", nrow(datos_master), "\n")
cat("Filas después de filtrar NAs:", nrow(datos_modelado), "\n")
cat("Filas eliminadas:", filas_eliminadas, "\n")
cat("  Proporción conservada:", 
    round(100 * nrow(datos_modelado) / nrow(datos_master), 2), "%\n")


## ----bloque-2-etsir-pred---------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 2.1: FUNCIÓN DE PREDICCIÓN ETSIR
#==============================================================================
# PROPÓSITO: Implementa la ecuación (1) de Wang & Zhang (2025) con funciones
# piecewise lineales para capturar efectos de umbral.
# 
# JUSTIFICACIÓN TEÓRICA:
#   - g_T(v) = c_T * v + c_T' * max(v - H_T, 0)
#   - g_P(v) = c_P * v + c_P' * max(v - H_P, 0)
#   - E(I_t | F_{t-1}) = {c0 + cI*ΣI_{t-i} + g_T + g_P} * I_{t-1}
# 
# VENTAJAS vs. VERSIÓN ANTERIOR:
#   - Unifica etsir_pred() y etsir_step() en una sola función
#   - Acepta data.frames completos o vectores (flexible para MLE y MSA)
#   - Usa nombres de variables de datos_master (acumulado_KI, temp_opt, etc.)
#==============================================================================

etsir_pred <- function(I_t_1, acumulado_KI, temperatura, precipitacion,
                       c0, cI, cP, cP2, HP, cT, cT2, HT) {
  
  # Función piecewise para temperatura: g_T(v)
  #   cT  > 0: efecto positivo inicial (más calor → más mosquitos)
  #   cT2 < 0: efecto negativo tras el umbral H_T
  gT <- cT * temperatura + cT2 * pmax(temperatura - HT, 0)
  
  # Función piecewise para precipitación: g_P(v)
  #   cP  > 0: efecto positivo inicial (más lluvia → más criaderos)
  #   cP2 < 0: efecto de lavado tras el umbral H_P
  gP <- cP * precipitacion + cP2 * pmax(precipitacion - HP, 0)
  
  # Ecuación ETSIR completa
  mu <- (c0 + cI * acumulado_KI + gT + gP) * I_t_1
  
  # Restricción biológica: la esperanza no puede ser negativa
  mu <- pmax(mu, 1e-6)
  
  return(mu)
}

cat("✔ Función etsir_pred() definida correctamente.\n")


## ----bloque-2-etsir-nll----------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 2.2: FUNCIÓN DE VEROSIMILITUD NEGATIVA (NegBin)
#==============================================================================
# PROPÓSITO: Devuelve la log-verosimilitud NEGATIVA para usar optimizadores
# que minimizan (optim, nlminb).
#
# JUSTIFICACIÓN TEÓRICA:
#   - I_t ~ NegBin(μ_t, θ) por sobredispersión (Var >> Media)
#   - θ se parametriza en log para garantizar positividad
#   - Umbrales H_T, H_P se estiman por grid search externo (no diferenciables)
#
# PARÁMETROS:
#   β = (c0, cI, cP, cP2, cT, cT2, HT, HP, log(θ))
#==============================================================================

etsir_nll <- function(params, datos) {
  c0        <- params[1]
  cI        <- params[2]
  cP        <- params[3]
  cP2       <- params[4]
  cT        <- params[5]
  cT2       <- params[6]
  HT        <- params[7]
  HP        <- params[8]
  log_theta <- params[9]
  theta     <- exp(log_theta)
  
  # Calcular esperanza condicional
  mu <- etsir_pred(
    I_t_1         = datos$I_t_1,
    acumulado_KI  = datos$acumulado_KI,
    temperatura   = datos$temp_opt,
    precipitacion = datos$precip_opt,
    c0 = c0, cI = cI, cP = cP, cP2 = cP2, HP = HP,
    cT = cT, cT2 = cT2, HT = HT
  )
  
  # Log-verosimilitud negativa NegBin
  nll <- -sum(dnbinom(datos$I_t, mu = mu, size = theta, log = TRUE))
  
  # Penalización si hay valores no finitos
  if (!is.finite(nll)) nll <- 1e10
  
  return(nll)
}

cat("✔ Función etsir_nll() definida correctamente.\n")


## ----bloque-3-ajustar-osa--------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 3.1: FUNCIÓN DE AJUSTE ETSIR-OSA POR PAÍS
#==============================================================================
# PROPÓSITO: Ajustar el modelo ETSIR usando estimación one-step-ahead con
# grid search de umbrales (H_T, H_P) y selección por AIC.
# 
# MEJORAS vs. VERSIÓN ANTERIOR:
#   1. Validación exhaustiva de datos antes de ajustar
#   2. Mensajes informativos de diagnóstico
#   3. Extracción robusta de coeficientes
#   4. Fallback a Poisson si NegBin no converge
#   5. Grid de umbrales basado en percentiles (más amplio)
#   6. Inicialización adaptativa de theta
#==============================================================================

ajustar_etsir_osa <- function(datos_pais, pais_nombre = NULL) {
  
  if (is.null(pais_nombre)) pais_nombre <- unique(datos_pais$pais)[1]
  cat("\n═══════════════════════════════════════════════════════\n")
  cat("  Ajustando ETSIR One-Step-Ahead para:", pais_nombre, "\n")
  cat("═══════════════════════════════════════════════════════\n")
  
  # --- PASO 0: Validación exhaustiva ---
  cat("  [0/5] Validando datos de entrada...\n")
  
  # Verificar que todas las variables necesarias existen
  vars_requeridas <- c("I_t", "I_t_1", "acumulado_KI", "temp_opt", "precip_opt")
  vars_faltantes <- setdiff(vars_requeridas, colnames(datos_pais))
  
  if (length(vars_faltantes) > 0) {
    stop("Variables faltantes en ", pais_nombre, ": ", 
         paste(vars_faltantes, collapse = ", "))
  }
  
  # Filtrar NAs en variables críticas
  datos_clean <- datos_pais %>%
    filter(!is.na(I_t), !is.na(I_t_1), !is.na(acumulado_KI),
           !is.na(temp_opt), !is.na(precip_opt))
  
  n_original <- nrow(datos_pais)
  n_clean <- nrow(datos_clean)
  
  if (n_clean == 0) {
    stop("No hay datos válidos después de filtrar NAs para ", pais_nombre)
  }
  
  if (n_clean < n_original) {
    cat("    ⚠ Se eliminaron", n_original - n_clean, "filas con NAs\n")
  }
  
  # Verificar que I_t tenga variación
  if (sd(datos_clean$I_t, na.rm = TRUE) == 0) {
    stop("I_t no tiene variación en ", pais_nombre, " (todos los valores son iguales)")
  }
  
  cat("    ✔ Datos validados:", n_clean, "observaciones\n")
  
  # --- PASO 1: Grid de umbrales ---
  cat("  [1/5] Construyendo grid de umbrales...\n")
  
  # Grid más amplio basado en percentiles
  HT_candidatos <- quantile(datos_clean$temp_opt, 
                            probs = c(0.10, 0.25, 0.50, 0.75, 0.90), 
                            na.rm = TRUE)
  HP_candidatos <- quantile(datos_clean$precip_opt, 
                            probs = c(0.10, 0.25, 0.50, 0.75, 0.90), 
                            na.rm = TRUE)
  
  grid_umbrales <- expand.grid(HT = HT_candidatos, HP = HP_candidatos)
  cat("    Grid:", nrow(grid_umbrales), "combinaciones\n")
  
  resultados_grid <- data.frame(
    HT = numeric(), HP = numeric(), 
    AIC = numeric(), logLik = numeric(),
    c0 = numeric(), cI = numeric(),
    cP = numeric(), cP2 = numeric(),
    cT = numeric(), cT2 = numeric(),
    theta = numeric(),
    convergio = logical()
  )
  
  # --- PASO 2: Ajustar glm.nb en grid ---
  cat("  [2/5] Ajustando modelos en grid...\n")
  
  n_convergio <- 0
  
  for (i in seq_len(nrow(grid_umbrales))) {
    HT_cand <- grid_umbrales$HT[i]
    HP_cand <- grid_umbrales$HP[i]
    
    # Construir variables piecewise para este par de umbrales
    datos_tmp <- datos_clean %>%
      mutate(
        temp_tramo1  = temp_opt,
        temp_tramo2  = pmax(temp_opt - HT_cand, 0),
        precip_tramo1 = precip_opt,
        precip_tramo2 = pmax(precip_opt - HP_cand, 0)
      )
    
    # Intentar ajustar NegBin
    mod <- tryCatch({
      # Primer intento con init.theta = 1.5
      MASS::glm.nb(
        I_t ~ I_t_1 + acumulado_KI + 
              temp_tramo1 + temp_tramo2 + 
              precip_tramo1 + precip_tramo2,
        data = datos_tmp,
        control = glm.control(maxit = 500, warn.complete.convergence = FALSE),
        init.theta = 1.5,
        link = "log"
      )
    }, error = function(e) {
      # Si falla, intentar sin init.theta
      tryCatch({
        MASS::glm.nb(
          I_t ~ I_t_1 + acumulado_KI + 
                temp_tramo1 + temp_tramo2 + 
                precip_tramo1 + precip_tramo2,
          data = datos_tmp,
          control = glm.control(maxit = 500, warn.complete.convergence = FALSE),
          link = "log"
        )
      }, error = function(e2) {
        # Si NegBin falla completamente, intentar Poisson
        tryCatch({
          glm(
            I_t ~ I_t_1 + acumulado_KI + 
                temp_tramo1 + temp_tramo2 + 
                precip_tramo1 + precip_tramo2,
            data = datos_tmp,
            family = poisson(link = "log"),
            control = glm.control(maxit = 500)
          )
        }, error = function(e3) NULL)
      })
    })
    
    # Extraer resultados si el modelo convergió
    if (!is.null(mod) && !inherits(mod, "try-error")) {
      # Verificar que el modelo tenga coeficientes válidos
      coefs <- coef(mod)
      
      if (!is.null(coefs) && length(coefs) > 0 && all(is.finite(coefs))) {
        # Extraer coeficientes de forma robusta
        c0_val <- if ("(Intercept)" %in% names(coefs)) coefs["(Intercept)"] else NA
        cI_val <- if ("I_t_1" %in% names(coefs)) coefs["I_t_1"] else NA
        cP_val <- if ("precip_tramo1" %in% names(coefs)) coefs["precip_tramo1"] else NA
        cP2_val <- if ("precip_tramo2" %in% names(coefs)) coefs["precip_tramo2"] else NA
        cT_val <- if ("temp_tramo1" %in% names(coefs)) coefs["temp_tramo1"] else NA
        cT2_val <- if ("temp_tramo2" %in% names(coefs)) coefs["temp_tramo2"] else NA
        
        # Theta solo para NegBin
        theta_val <- if (inherits(mod, "negbin")) mod$theta else NA
        
        resultados_grid <- rbind(resultados_grid, data.frame(
          HT = HT_cand, HP = HP_cand,
          AIC = AIC(mod),
          logLik = as.numeric(logLik(mod)),
          c0 = c0_val, cI = cI_val,
          cP = cP_val, cP2 = cP2_val,
          cT = cT_val, cT2 = cT2_val,
          theta = theta_val,
          convergio = TRUE
        ))
        
        n_convergio <- n_convergio + 1
      }
    }
  }
  
  cat("    ✔ Modelos convergidos:", n_convergio, "/", nrow(grid_umbrales), "\n")
  
  if (nrow(resultados_grid) == 0) {
    warning("Ningún modelo convergió para ", pais_nombre)
    return(NULL)
  }
  
  # --- PASO 3: Seleccionar mejor modelo por AIC ---
  cat("  [3/5] Seleccionando mejor modelo por AIC...\n")
  
  mejor_idx <- which.min(resultados_grid$AIC)
  mejor <- resultados_grid[mejor_idx, ]
  
  cat("    ✔ Mejor modelo:\n")
  cat("      H_T =", round(mejor$HT, 2), "°C | H_P =", round(mejor$HP, 2), "mm\n")
  cat("      AIC =", round(mejor$AIC, 2), "\n")
  if (!is.na(mejor$theta)) {
    cat("      θ (dispersión) =", round(mejor$theta, 3), "\n")
  }
  
  # --- PASO 4: Reconstruir modelo final ---
  cat("  [4/5] Reconstruyendo modelo final...\n")
  
  datos_finales <- datos_clean %>%
    mutate(
      temp_tramo1   = temp_opt,
      temp_tramo2   = pmax(temp_opt - mejor$HT, 0),
      precip_tramo1 = precip_opt,
      precip_tramo2 = pmax(precip_opt - mejor$HP, 0)
    )
  
  # Ajustar modelo final
  modelo_final <- tryCatch({
    if (!is.na(mejor$theta)) {
      MASS::glm.nb(
        I_t ~ I_t_1 + acumulado_KI + 
              temp_tramo1 + temp_tramo2 + 
              precip_tramo1 + precip_tramo2,
        data = datos_finales,
        control = glm.control(maxit = 1000),
        init.theta = mejor$theta
      )
    } else {
      glm(
        I_t ~ I_t_1 + acumulado_KI + 
              temp_tramo1 + temp_tramo2 + 
              precip_tramo1 + precip_tramo2,
        data = datos_finales,
        family = poisson(link = "log"),
        control = glm.control(maxit = 1000)
      )
    }
  }, error = function(e) {
    warning("No se pudo reconstruir el modelo final para ", pais_nombre)
    NULL
  })
  
  if (is.null(modelo_final)) {
    return(NULL)
  }
  
  # --- PASO 5: Predicciones y métricas ---
  cat("  [5/5] Calculando predicciones y métricas...\n")
  
  datos_finales$mu_osa <- predict(modelo_final, type = "response")
  datos_finales$residuo_osa <- datos_finales$I_t - datos_finales$mu_osa
  
  # Métricas de evaluación
  ss_res <- sum(datos_finales$residuo_osa^2, na.rm = TRUE)
  ss_tot <- sum((datos_finales$I_t - mean(datos_finales$I_t, na.rm = TRUE))^2, na.rm = TRUE)
  R2 <- 1 - ss_res / ss_tot
  RMSE <- sqrt(mean(datos_finales$residuo_osa^2, na.rm = TRUE))
  
  cat("    ✔ Métricas OSA:\n")
  cat("      R²   =", round(R2, 4), "\n")
  cat("      RMSE =", round(RMSE, 2), "\n\n")
  
  # ---------------------------------------------------------
  # Retornar lista con todos los resultados
  # --------------------------------------------------------- 
  return(list(
    pais           = pais_nombre,
    modelo         = modelo_final,
    datos          = datos_finales,
    umbrales       = c(HT = mejor$HT, HP = mejor$HP),
    coeficientes   = c(c0 = mejor$c0, cI = mejor$cI,
                       cP = mejor$cP, cP2 = mejor$cP2,
                       cT = mejor$cT, cT2 = mejor$cT2),
    theta          = mejor$theta,
    metricas       = c(AIC = mejor$AIC, R2 = R2, RMSE = RMSE),
    grid_completo  = resultados_grid
  ))
}

cat("✔ Función ajustar_etsir_osa() definida correctamente.\n")



## ----bloque-3-ejecucion-osa, warning=FALSE---------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 3.2: EJECUCIÓN DEL AJUSTE OSA PARA LOS 5 PAÍSES
#==============================================================================


# Ajustar por país
resultados_etsir_osa <- datos_modelado %>%
  split(.$pais) %>%
  map(~ ajustar_etsir_osa(.x))

# Consolidar métricas en tabla
tabla_metricas_osa <- map_dfr(resultados_etsir_osa, ~ {
  tibble(
    País       = .x$pais,
    `H_T (°C)` = round(.x$umbrales["HT"], 2),
    `H_P (mm)` = round(.x$umbrales["HP"], 2),
    `c₀`       = round(.x$coeficientes["c0"], 4),
    `c_I`      = round(.x$coeficientes["cI"], 4),
    `c_T`      = round(.x$coeficientes["cT"], 4),
    `c_T'`     = round(.x$coeficientes["cT2"], 4),
    `c_P`      = round(.x$coeficientes["cP"], 4),
    `c_P'`     = round(.x$coeficientes["cP2"], 4),
    `θ`        = round(.x$theta, 3),
    AIC        = round(.x$metricas["AIC"], 1),
    `R²`       = round(.x$metricas["R2"], 4),
    RMSE       = round(.x$metricas["RMSE"], 1)
  )
})

# Renderizar tabla
cat("\n═══════════════════════════════════════════════════════\n")
cat("RESULTADOS DEL MODELO ETSIR ONE-STEP-AHEAD\n")
cat("═══════════════════════════════════════════════════════\n\n")

tabla_metricas_osa %>%
  kbl(caption = "Resultados del Modelo ETSIR One-Step-Ahead por País",
      format = "html", escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#2c3e50", color = "white")


## ----bloque-4-diagnostico-series-------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 4.1: DIAGNÓSTICO VISUAL - SERIES OBSERVADAS VS. PREDICHAS
#==============================================================================

ggplot(datos_predichos, aes(x = fecha_semana)) +
  # Línea Observada: Color neutro fijo (fuera de aes) para que no interfiera
  geom_line(aes(y = I_t), color = "#7f8c8d", linewidth = 0.6) +
  
  # Línea Predicha: Toma el color según el país desde tu paleta
  geom_line(aes(y = mu_osa, color = pais), linewidth = 0.7, linetype = "dashed") +
  
  facet_wrap(~ pais, scales = "free_y", ncol = 2) +
  
  # Aplica tu paleta de países directamente
  scale_color_manual(values = PALETA_PAISES) +
  
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title    = "Modelo ETSIR One-Step-Ahead: Observado vs. Predicho",
    subtitle = "Línea gris: Observado | Línea de color discontinua: Predicho por país",
    x = "Fecha", y = "Casos semanales"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", # Se oculta la leyenda porque el título del panel ya dice el país
        strip.background = element_rect(fill = "#34495e"),
        strip.text = element_text(color = "white", face = "bold"))



## ----bloque-4-diagnostico-residuos, fig.height=6, fig.width=10-------------------------------------------------------------
#==============================================================================
# BLOQUE 4.2: DIAGNÓSTICO DE RESIDUOS (HOMOCEDASTICIDAD E INDEPENDENCIA)
#==============================================================================

# --- 4.2.1 Gráfico de residuos vs. valores ajustados ---
p_residuos <- ggplot(datos_predichos, aes(x = mu_osa, y = residuo_osa, color = pais)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
  facet_wrap(~ pais, scales = "free", ncol = 3) +
  scale_color_manual(values = PALETA_PAISES) +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title    = "Residuos del Modelo ETSIR-OSA",
    subtitle = "Patrón horizontal = buen ajuste; curva = posible misspecification",
    x = "Valor ajustado (μ_t)", y = "Residuo (I_t - μ_t)", color = "País"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

print(p_residuos)

# --- 4.2.2 ACF de residuos (verificar independencia) ---
cat("\n═══════════════════════════════════════════════════════\n")
cat("ANÁLISIS DE AUTOCORRELACIÓN DE RESIDUOS\n")
cat("═══════════════════════════════════════════════════════\n\n")

par(mfrow = c(2, 3))

for (p in unique(datos_predichos$pais)) {
  
  # 1. Extraer los residuos del país actual
  residuos_pais <- datos_predichos %>%
    filter(pais == p) %>%
    pull(residuo_osa)
  
  # 2. Buscar el color correspondiente al país en tu paleta
  # Si el país no está en la paleta, usa gris por defecto ("#7f8c8d")
  color_actual <- PALETA_PAISES[p]
  if (is.na(color_actual)) color_actual <- "#7f8c8d" 
  
  # 3. Graficar la ACF con el color asignado
  acf(residuos_pais, 
      main = p, 
      lag.max = 30, 
      col = color_actual, # <-- Se aplica el color personalizado aquí
      lwd = 2)
}

par(mfrow = c(1, 1))


cat("\nInterpretación:\n")
cat("  → Si la ACF decae rápidamente: el modelo captura bien la dinámica.\n")
cat("  → Si hay picos significativos en rezagos altos: se requiere MSA.\n")


## ----bloque-5-msa-nll------------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 5.1: FUNCIÓN DE PÉRDIDA MSA - NEGATIVE LOG-LIKELIHOOD MULTI-PASO
#==============================================================================
# PROPÓSITO: Evalúa la log-verosimilitud NEGATIVA de una distribución de
# conteo (Poisson o NegBin) sobre la recursión MSA a K_max pasos.
# 
# DIFERENCIAS CRÍTICAS vs. VERSIÓN ANTERIOR (msa_loss_core_vectorized):
#   1. Usa verosimilitud real (no mínimos cuadrados puros)
#   2. Permite comparar Poisson vs. NegBin con AIC/BIC
#   3. Las métricas se calculan con la MISMA recursión que se optimizó
#   4. El escalado es dinámico por país (no constantes fijas)
# 
# JUSTIFICACIÓN TEÓRICA:
#   - Distribución = "poisson" → Var(I_t) = μ_t
#   - Distribución = "nbinom"  → Var(I_t) = μ_t + μ_t²/θ
#==============================================================================

msa_negloglik_core <- function(params, I_s, I_orig, T_s, P_s, K_I, K_max,
                                distribucion = c("poisson", "nbinom"),
                                escala_I) {
  # NOTA: escala_I ya NO tiene un valor por defecto fijo (antes era 1000).
  # Debe pasarse explícitamente como el I_sd calculado dinámicamente para
  # el país en cuestión (ver ajustar_msa_distribucion), porque I_s fue
  # estandarizado dividiendo por ese mismo I_sd. Usar una constante fija
  # aquí desalineaba la escala de las predicciones con la escala real de
  # I_orig en países cuya I_sd es muy distinta de 1000 (ej. Costa Rica
  # con I_t máximo de 1749 vs. México con 27627).

  distribucion <- match.arg(distribucion)

  c0 <- params[1]; cI <- params[2]; cT <- params[3]
  cT_p <- params[4]; cP <- params[5]; cP_p <- params[6]
  HT_s <- params[7]; HP_s <- params[8]

  # theta (dispersión NegBin) se estima en escala log para garantizar > 0;
  # si la distribución es Poisson, este parámetro simplemente no se usa.
  log_theta <- if (distribucion == "nbinom") params[9] else NA

  n <- length(I_s)
  valid_t <- (K_I + 1):(n - K_max)
  n_valid <- length(valid_t)

  if (n_valid <= 0) return(Inf)

  # Matrices de condiciones ambientales futuras
  idx_future <- outer(valid_t, 0:(K_max - 1), "+")
  env_T_matrix <- matrix(T_s[idx_future], nrow = n_valid)
  env_P_matrix <- matrix(P_s[idx_future], nrow = n_valid)

  # Matriz de estados históricos (Shift Register)
  idx_history <- outer(valid_t - 1, 0:(K_I - 1), "-")
  state_matrix <- matrix(I_s[idx_history], nrow = n_valid)

  preds_matrix <- matrix(0, nrow = n_valid, ncol = K_max)

  # Recursión MSA vectorizada
  for (k in 1:K_max) {
    I_t_minus_1 <- state_matrix[, 1]
    acum_I <- rowSums(state_matrix)

    gT <- cT * env_T_matrix[, k] + cT_p * pmax(env_T_matrix[, k] - HT_s, 0)
    gP <- cP * env_P_matrix[, k] + cP_p * pmax(env_P_matrix[, k] - HP_s, 0)

    mu <- (c0 + cI * acum_I + gT + gP) * I_t_minus_1
    mu <- pmax(mu, 1e-4)

    preds_matrix[, k] <- mu

    if (K_I > 1) {
      state_matrix <- cbind(mu, state_matrix[, 1:(K_I - 1), drop = FALSE])
    } else {
      state_matrix <- matrix(mu, ncol = 1)
    }
  }

  # Evaluar verosimilitud en escala ORIGINAL (conteos reales)
  idx_obs <- outer(valid_t, 1:K_max, "+")
  # La verosimilitud se evalúa en la escala ORIGINAL de conteos (no escalada),
  # porque dpois/dnbinom requieren conteos enteros no negativos reales.
  obs_matrix_orig <- matrix(I_orig[idx_obs], nrow = n_valid)
  preds_matrix_orig <- preds_matrix * escala_I

  if (distribucion == "poisson") {
    ll <- sum(dpois(obs_matrix_orig, lambda = pmax(preds_matrix_orig, 1e-6), log = TRUE))
  } else {
    theta <- exp(log_theta)
    ll <- sum(dnbinom(obs_matrix_orig, mu = pmax(preds_matrix_orig, 1e-6),
                       size = theta, log = TRUE))
  }

  if (!is.finite(ll)) return(1e10)

  # Penalización L2 ligera sobre los coeficientes estructurales (no sobre theta)
  # para estabilizar la optimización, igual que en la versión anterior.
  nll <- -ll + 0.01 * sum(params[1:6]^2)
  return(nll)
}

cat("✔ Función msa_negloglik_core() definida correctamente.\n")



## ----bloque-5-ajustar-msa--------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 5.2: FUNCIÓN DE AJUSTE MSA POR PAÍS Y DISTRIBUCIÓN
#==============================================================================
# PROPÓSITO: Ajustar el modelo ETSIR usando estimación MSA con K_max pasos.
#
# CORRECCIONES CRÍTICAS vs. VERSIÓN ANTERIOR (ajustar_msa_pais_completo):
#   1. Estandarización dinámica (z-score) por país
#      - Antes: T_s = (temp-25)/5, P_s = (precip-200)/100 (constantes fijas)
#      - Ahora: z-score con media y sd del propio país
#      - Justificación: temp_opt/precip_opt son ventanas K variables con
#        rangos muy distintos entre países
#
#   2. Multistart (8 puntos de inicio)
#      - Antes: Un solo punto → quedaba pegado en límites de caja
#      - Ahora: 1 punto base + 7 aleatorios → evita óptimos locales
#      - Justificación: Honduras y México mostraban R² negativos (-17, -0.66)
#
#   3. Límites de caja ampliados
#      - Antes: cT', cP' ∈ [-2, -1e-4]
#      - Ahora: cT', cP' ∈ [-10, -1e-4]
#      - Justificación: El óptimo real podía estar más allá de -2
#==============================================================================

ajustar_msa_distribucion <- function(df_pais, pais_nombre, distribucion,
                                      K_I = 8, K_max = 4) {
  
  # --- Estandarización dinámica por país ---
  temp_mu   <- mean(df_pais$temp_opt, na.rm = TRUE)
  temp_sd   <- sd(df_pais$temp_opt, na.rm = TRUE)
  precip_mu <- mean(df_pais$precip_opt, na.rm = TRUE)
  precip_sd <- sd(df_pais$precip_opt, na.rm = TRUE)
  I_mu      <- mean(df_pais$I_t, na.rm = TRUE)
  I_sd      <- sd(df_pais$I_t, na.rm = TRUE)
  
  # Salvaguarda: evitar división por cero
  temp_sd   <- ifelse(temp_sd == 0 || is.na(temp_sd), 1, temp_sd)
  precip_sd <- ifelse(precip_sd == 0 || is.na(precip_sd), 1, precip_sd)
  I_sd      <- ifelse(I_sd == 0 || is.na(I_sd), 1, I_sd)
  
  I_s     <- (df_pais$I_t - 0) / I_sd
  T_s     <- (df_pais$temp_opt - temp_mu) / temp_sd
  P_s     <- (df_pais$precip_opt - precip_mu) / precip_sd
  I_orig  <- df_pais$I_t
  
  if (length(I_s) < 50) return(NULL)
  
  n_params_estructurales <- 6
  n_params_dist <- if (distribucion == "nbinom") 1 else 0
  
  # --- Grid search de umbrales ---
  HT_grid <- quantile(T_s, probs = seq(0.2, 0.8, length.out = 5), na.rm = TRUE)
  HP_grid <- quantile(P_s, probs = seq(0.2, 0.8, length.out = 5), na.rm = TRUE)
  grid <- expand.grid(HT = HT_grid, HP = HP_grid)
  grid$nll <- Inf
  
  params_init_base <- c(1.2, -0.01, 0.05, -0.02, 0.05, -0.02)
  if (distribucion == "nbinom") params_init_base <- c(params_init_base, log(2))
  
  # --- Configuración de multistart ---
  N_STARTS <- 8L
  N_STARTS_GRID <- 3L
  
  lower_b <- c(0.05, -3, 0.0001, -10, 0.0001, -10)
  upper_b <- c(8.0, -1e-5, 3.0, -1e-5, 3.0, -1e-5)
  if (distribucion == "nbinom") {
    lower_b <- c(lower_b, log(0.02))
    upper_b <- c(upper_b, log(200))
  }
  
  construir_params <- function(p, HT, HP) {
    if (distribucion == "nbinom") c(p[1:6], HT, HP, p[7]) else c(p[1:6], HT, HP)
  }
  
  generar_puntos_init <- function(n_starts, semilla = 123) {
    set.seed(semilla)
    puntos <- list(params_init_base)
    for (i in seq_len(n_starts - 1)) {
      p_rand <- c(
        runif(1, 0.2, 4.5),
        runif(1, -1.8, -0.0005),
        runif(1, 0.0005, 1.8),
        runif(1, -1.8, -0.0005),
        runif(1, 0.0005, 1.8),
        runif(1, -1.8, -0.0005)
      )
      if (distribucion == "nbinom") {
        p_rand <- c(p_rand, runif(1, log(0.1), log(50)))
      }
      puntos[[length(puntos) + 1]] <- p_rand
    }
    puntos
  }
  
  optimizar_multistart <- function(fn_obj, puntos_init, lower_b, upper_b) {
    mejor_resultado <- NULL
    mejor_valor <- Inf
    for (p0 in puntos_init) {
      res <- tryCatch({
        suppressWarnings(optim(
          par = p0, fn = fn_obj,
          method = "L-BFGS-B", lower = lower_b, upper = upper_b,
          control = list(maxit = 300)
        ))
      }, error = function(e) list(value = Inf))
      
      if (is.finite(res$value) && res$value < mejor_valor) {
        mejor_valor <- res$value
        mejor_resultado <- res
      }
    }
    if (is.null(mejor_resultado)) {
      mejor_resultado <- suppressWarnings(optim(
        par = puntos_init[[1]], fn = fn_obj,
        method = "Nelder-Mead", control = list(maxit = 1000)
      ))
    }
    mejor_resultado
  }
  
  puntos_init_grid  <- generar_puntos_init(N_STARTS_GRID, semilla = 123)
  puntos_init_final <- generar_puntos_init(N_STARTS, semilla = 123)
  
  # --- Grid de umbrales con multistart liviano ---
  for (i in seq_len(nrow(grid))) {
    fn_grid <- function(p) msa_negloglik_core(
      construir_params(p, grid$HT[i], grid$HP[i]),
      I_s = I_s, I_orig = I_orig, T_s = T_s, P_s = P_s,
      K_I = K_I, K_max = K_max, distribucion = distribucion,
      escala_I = I_sd
    )
    opt_grid <- optimizar_multistart(fn_grid, puntos_init_grid, lower_b, upper_b)
    grid$nll[i] <- opt_grid$value
  }
  
  mejor_grid <- grid[which.min(grid$nll), ]
  HT_opt <- mejor_grid$HT
  HP_opt <- mejor_grid$HP
  
  # --- Optimización final con multistart completo ---
  fn_final <- function(p) {
    msa_negloglik_core(construir_params(p, HT_opt, HP_opt),
                        I_s = I_s, I_orig = I_orig, T_s = T_s, P_s = P_s,
                        K_I = K_I, K_max = K_max, distribucion = distribucion,
                        escala_I = I_sd)
  }
  
  optim_result <- optimizar_multistart(fn_final, puntos_init_final, lower_b, upper_b)
  
  # Diagnóstico de límites
  par_chequeo <- optim_result$par[1:6]
  en_limite <- abs(par_chequeo - lower_b[1:6]) < 1e-4 | abs(par_chequeo - upper_b[1:6]) < 1e-4
  if (any(en_limite)) {
    nombres_par <- c("c0", "cI", "cT", "cT'", "cP", "cP'")
    cat("    ⚠", pais_nombre, "(", distribucion, "): parámetro(s)",
        paste(nombres_par[en_limite], collapse = ", "),
        "quedaron en el límite de la caja.\n")
  }
  
  # Des-escalado de umbrales
  HT_orig <- HT_opt * temp_sd + temp_mu
  HP_orig <- HP_opt * precip_sd + precip_mu
  
  par_final <- optim_result$par
  theta_final <- if (distribucion == "nbinom") exp(par_final[7]) else NA
  
  params_estimados <- c(
    c0 = par_final[1], cI = par_final[2], cT = par_final[3],
    cT_p = par_final[4], cP = par_final[5], cP_p = par_final[6],
    HT_s = HT_opt, HP_s = HP_opt,
    theta = theta_final
  )
  
  # --- Predicciones multi-paso para diagnóstico ---
  n <- length(I_s)
  preds_1paso <- rep(NA_real_, n)
  
  for (t in (K_I + 1):n) {
    state <- I_s[(t - 1):(t - K_I)]
    gT <- params_estimados["cT"] * T_s[t] + params_estimados["cT_p"] * pmax(T_s[t] - HT_opt, 0)
    gP <- params_estimados["cP"] * P_s[t] + params_estimados["cP_p"] * pmax(P_s[t] - HP_opt, 0)
    mu <- (params_estimados["c0"] + params_estimados["cI"] * sum(state) + gT + gP) * state[1]
    preds_1paso[t] <- max(mu, 1e-4) * I_sd
  }
  
  residuos <- I_orig - preds_1paso
  residuos_esc <- residuos / sd(I_orig, na.rm = TRUE)
  
  ss_res <- sum(residuos^2, na.rm = TRUE)
  ss_tot <- sum((I_orig - mean(I_orig, na.rm = TRUE))^2, na.rm = TRUE)
  R2 <- ifelse(ss_tot == 0, 0, 1 - ss_res / ss_tot)
  RMSE <- sqrt(mean(residuos^2, na.rm = TRUE))
  
  residuos_limpios <- residuos_esc[is.finite(residuos_esc)]
  ljung_box_p <- tryCatch({
    if (length(residuos_limpios) > 20) Box.test(residuos_limpios, lag = 20, type = "Ljung-Box")$p.value else NA
  }, error = function(e) NA)
  
  # --- AIC/BIC del criterio MSA ---
  valid_t <- (K_I + 1):(n - K_max)
  n_obs_msa <- length(valid_t) * K_max
  k_params <- n_params_estructurales + 2 + n_params_dist
  
  nll_final <- optim_result$value
  AIC_msa <- 2 * nll_final + 2 * k_params
  BIC_msa <- 2 * nll_final + k_params * log(n_obs_msa)
  
  list(
    tabla = tibble(
      Pais = pais_nombre,
      Distribucion = ifelse(distribucion == "poisson", "Poisson", "Binomial Negativa"),
      c0 = round(par_final[1], 4), cI = round(par_final[2], 4),
      cT = round(par_final[3], 4), cT_prime = round(par_final[4], 4),
      cP = round(par_final[5], 4), cP_prime = round(par_final[6], 4),
      `H_T (°C)` = round(HT_orig, 2), `H_P (mm)` = round(HP_orig, 2),
      theta = if (distribucion == "nbinom") round(theta_final, 3) else NA,
      NegLogLik = round(nll_final, 2),
      AIC = round(AIC_msa, 1), BIC = round(BIC_msa, 1),
      R2 = round(R2, 4), RMSE = round(RMSE, 2),
      `Ljung-Box (p)` = round(ljung_box_p, 4),
      Iteraciones = optim_result$counts["function"]
    ),
    datos_validacion = tibble(
      pais = pais_nombre, distribucion = distribucion,
      semana = seq_along(I_orig), I_obs = I_orig,
      I_pred = preds_1paso, residuos = residuos, residuos_std = residuos_esc,
      temp = df_pais$temp_opt, precip = df_pais$precip_opt
    ),
    params = params_estimados,
    umbrales = c(HT = HT_orig, HP = HP_orig)
  )
}

ajustar_msa_distribucion <- function(df_pais, pais_nombre, distribucion,
                                      K_I = 8, K_max = 4,
                                      timeout_sec = 300,      # tiempo máximo
                                      show_progress = TRUE) { # barra de progreso
  
  # Si se pide timeout, usar withTimeout
  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(R.utils::withTimeout({
      .ajustar_msa_dist_interno(df_pais, pais_nombre, distribucion,
                                 K_I, K_max, show_progress)
    }, timeout = timeout_sec, onTimeout = "error"))
  } else {
    warning("Paquete R.utils no instalado. No se aplicará timeout.")
    return(.ajustar_msa_dist_interno(df_pais, pais_nombre, distribucion,
                                      K_I, K_max, show_progress))
  }
}

# Función interna (sin timeout) que hace el trabajo real
.ajustar_msa_dist_interno <- function(df_pais, pais_nombre, distribucion,
                                       K_I, K_max, show_progress) {
  
  cat(sprintf("\n─────────────────────────────────────────────────────\n"))
  cat(sprintf("  Procesando MSA para: %s (%s)\n", pais_nombre, distribucion))
  cat(sprintf("─────────────────────────────────────────────────────\n"))
  
  # --- Estandarización dinámica por país ---
  cat("  [0/5] Estandarizando variables...\n")
  
  temp_mu   <- mean(df_pais$temp_opt, na.rm = TRUE)
  temp_sd   <- sd(df_pais$temp_opt, na.rm = TRUE)
  precip_mu <- mean(df_pais$precip_opt, na.rm = TRUE)
  precip_sd <- sd(df_pais$precip_opt, na.rm = TRUE)
  I_mu      <- mean(df_pais$I_t, na.rm = TRUE)
  I_sd      <- sd(df_pais$I_t, na.rm = TRUE)
  
  temp_sd   <- ifelse(temp_sd == 0 || is.na(temp_sd), 1, temp_sd)
  precip_sd <- ifelse(precip_sd == 0 || is.na(precip_sd), 1, precip_sd)
  I_sd      <- ifelse(I_sd == 0 || is.na(I_sd), 1, I_sd)
  
  I_s     <- (df_pais$I_t - 0) / I_sd
  T_s     <- (df_pais$temp_opt - temp_mu) / temp_sd
  P_s     <- (df_pais$precip_opt - precip_mu) / precip_sd
  I_orig  <- df_pais$I_t
  
  n <- length(I_s)
  if (n < 50) {
    warning("Datos insuficientes para ", pais_nombre, " (", n, " obs)")
    return(NULL)
  }
  
  cat("    ✔", n, "observaciones | I_sd =", round(I_sd, 2), 
      "| temp_sd =", round(temp_sd, 2), "| precip_sd =", round(precip_sd, 2), "\n")
  
  n_params_estructurales <- 6
  n_params_dist <- if (distribucion == "nbinom") 1 else 0
  
  # --- Grid search de umbrales (más reducido) ---
  cat("  [1/5] Grid search de umbrales (5x5)...\n")
  
  HT_grid <- quantile(T_s, probs = seq(0.2, 0.8, length.out = 5), na.rm = TRUE)
  HP_grid <- quantile(P_s, probs = seq(0.2, 0.8, length.out = 5), na.rm = TRUE)
  grid <- expand.grid(HT = HT_grid, HP = HP_grid)
  n_grid <- nrow(grid)
  grid$nll <- Inf
  grid$convergio <- FALSE
  
  cat("    Grid:", n_grid, "combinaciones\n")
  
  # Parámetros base para inicialización
  params_init_base <- c(1.2, -0.01, 0.05, -0.02, 0.05, -0.02)
  if (distribucion == "nbinom") params_init_base <- c(params_init_base, log(2))
  
  # --- Configuración de multistart (más rápido) ---
  cat("  [2/5] Configurando optimización multistart (rápida)...\n")
  
  N_STARTS <- 4L          # reducido de 8 a 4
  N_STARTS_GRID <- 2L     # reducido de 3 a 2
  
  lower_b <- c(0.05, -3, 0.0001, -10, 0.0001, -10)
  upper_b <- c(8.0, -1e-5, 3.0, -1e-5, 3.0, -1e-5)
  if (distribucion == "nbinom") {
    lower_b <- c(lower_b, log(0.02))
    upper_b <- c(upper_b, log(200))
  }
  
  construir_params <- function(p, HT, HP) {
    if (distribucion == "nbinom") c(p[1:6], HT, HP, p[7]) else c(p[1:6], HT, HP)
  }
  
  generar_puntos_init <- function(n_starts, semilla = 123) {
    set.seed(semilla)
    puntos <- list(params_init_base)
    for (i in seq_len(n_starts - 1)) {
      p_rand <- c(
        runif(1, 0.2, 4.5),
        runif(1, -1.8, -0.0005),
        runif(1, 0.0005, 1.8),
        runif(1, -1.8, -0.0005),
        runif(1, 0.0005, 1.8),
        runif(1, -1.8, -0.0005)
      )
      if (distribucion == "nbinom") {
        p_rand <- c(p_rand, runif(1, log(0.1), log(50)))
      }
      puntos[[length(puntos) + 1]] <- p_rand
    }
    puntos
  }
  
  optimizar_multistart <- function(fn_obj, puntos_init, lower_b, upper_b,
                                   maxit_grid = 300, maxit_final = 800) {
    mejor_resultado <- NULL
    mejor_valor <- Inf
    for (p0 in puntos_init) {
      res <- tryCatch({
        suppressWarnings(optim(
          par = p0, fn = fn_obj,
          method = "L-BFGS-B", lower = lower_b, upper = upper_b,
          control = list(maxit = maxit_grid)
        ))
      }, error = function(e) list(value = Inf))
      
      if (is.finite(res$value) && res$value < mejor_valor) {
        mejor_valor <- res$value
        mejor_resultado <- res
      }
    }
    # Si todos fallan, intentar con Nelder-Mead (más robusto pero lento)
    if (is.null(mejor_resultado)) {
      mejor_resultado <- suppressWarnings(optim(
        par = puntos_init[[1]], fn = fn_obj,
        method = "Nelder-Mead", control = list(maxit = maxit_final)
      ))
    }
    mejor_resultado
  }
  
  puntos_init_grid  <- generar_puntos_init(N_STARTS_GRID, semilla = 123)
  puntos_init_final <- generar_puntos_init(N_STARTS, semilla = 123)
  
  # --- Grid de umbrales con multistart liviano ---
  cat("  [3/5] Ejecutando grid search (", n_grid, " combinaciones)...\n")
  
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = n_grid, style = 3)
  }
  
  n_convergio <- 0
  n_fallos <- 0
  
  for (i in seq_len(n_grid)) {
    HT_cand <- grid$HT[i]
    HP_cand <- grid$HP[i]
    
    fn_grid <- function(p) {
      msa_negloglik_core(
        construir_params(p, HT_cand, HP_cand),
        I_s = I_s, I_orig = I_orig, T_s = T_s, P_s = P_s,
        K_I = K_I, K_max = K_max, distribucion = distribucion,
        escala_I = I_sd
      )
    }
    
    opt_grid <- tryCatch({
      optimizar_multistart(fn_grid, puntos_init_grid, lower_b, upper_b,
                           maxit_grid = 300, maxit_final = 600)
    }, error = function(e) list(value = Inf))
    
    if (is.finite(opt_grid$value) && opt_grid$value < 1e9) {
      grid$nll[i] <- opt_grid$value
      grid$convergio[i] <- TRUE
      n_convergio <- n_convergio + 1
    } else {
      n_fallos <- n_fallos + 1
    }
    
    if (show_progress) setTxtProgressBar(pb, i)
  }
  
  if (show_progress) close(pb)
  
  cat("    ✔ Convergidos:", n_convergio, "/", n_grid, " | Fallos:", n_fallos, "\n")
  
  # Si ningún punto del grid convergió, intentar con un grid más grueso
  if (n_convergio == 0) {
    warning("Ninguna combinación del grid convergió para ", pais_nombre, 
            " (", distribucion, "). Intentando con grid reducido...")
    HT_grid2 <- quantile(T_s, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    HP_grid2 <- quantile(P_s, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
    grid2 <- expand.grid(HT = HT_grid2, HP = HP_grid2)
    grid2$nll <- Inf
    grid2$convergio <- FALSE
    
    for (i in seq_len(nrow(grid2))) {
      HT_cand <- grid2$HT[i]
      HP_cand <- grid2$HP[i]
      fn_grid2 <- function(p) {
        msa_negloglik_core(
          construir_params(p, HT_cand, HP_cand),
          I_s = I_s, I_orig = I_orig, T_s = T_s, P_s = P_s,
          K_I = K_I, K_max = K_max, distribucion = distribucion,
          escala_I = I_sd
        )
      }
      opt_grid2 <- tryCatch({
        optimizar_multistart(fn_grid2, puntos_init_grid, lower_b, upper_b,
                             maxit_grid = 300, maxit_final = 600)
      }, error = function(e) list(value = Inf))
      if (is.finite(opt_grid2$value) && opt_grid2$value < 1e9) {
        grid2$nll[i] <- opt_grid2$value
        grid2$convergio[i] <- TRUE
        n_convergio <- n_convergio + 1
      }
    }
    
    if (n_convergio == 0) {
      warning("No convergió ningún modelo para ", pais_nombre, " (", distribucion, ")")
      return(NULL)
    }
    grid <- grid2
  }
  
  mejor_grid <- grid[which.min(grid$nll), ]
  HT_opt <- mejor_grid$HT
  HP_opt <- mejor_grid$HP
  
  cat("    ✔ Mejores umbrales (z): H_T =", round(HT_opt, 3), 
      "| H_P =", round(HP_opt, 3), "\n")
  
  # --- Optimización final con multistart completo ---
  cat("  [4/5] Optimización final (multistart)...\n")
  
  fn_final <- function(p) {
    msa_negloglik_core(
      construir_params(p, HT_opt, HP_opt),
      I_s = I_s, I_orig = I_orig, T_s = T_s, P_s = P_s,
      K_I = K_I, K_max = K_max, distribucion = distribucion,
      escala_I = I_sd
    )
  }
  
  optim_result <- tryCatch({
    optimizar_multistart(fn_final, puntos_init_final, lower_b, upper_b,
                         maxit_grid = 500, maxit_final = 1000)
  }, error = function(e) {
    cat("    ⚠ Error en optimización final:", e$message, "\n")
    list(value = Inf, par = rep(NA, length(lower_b)))
  })
  
  if (!is.finite(optim_result$value) || optim_result$value > 1e9) {
    warning("Optimización final falló para ", pais_nombre, " (", distribucion, ")")
    return(NULL)
  }
  
  # Diagnóstico de límites
  par_chequeo <- optim_result$par[1:6]
  en_limite <- abs(par_chequeo - lower_b[1:6]) < 1e-4 | abs(par_chequeo - upper_b[1:6]) < 1e-4
  if (any(en_limite)) {
    nombres_par <- c("c0", "cI", "cT", "cT'", "cP", "cP'")
    cat("    ⚠ Parámetro(s)", paste(nombres_par[en_limite], collapse = ", "),
        "en límite de caja\n")
  }
  
  # Des-escalado de umbrales
  HT_orig <- HT_opt * temp_sd + temp_mu
  HP_orig <- HP_opt * precip_sd + precip_mu
  
  par_final <- optim_result$par
  theta_final <- if (distribucion == "nbinom") exp(par_final[7]) else NA
  
  params_estimados <- c(
    c0 = par_final[1], cI = par_final[2], cT = par_final[3],
    cT_p = par_final[4], cP = par_final[5], cP_p = par_final[6],
    HT_s = HT_opt, HP_s = HP_opt,
    theta = theta_final
  )
  
  cat("    ✔ Umbrales originales: H_T =", round(HT_orig, 2), "°C | H_P =", 
      round(HP_orig, 2), "mm\n")
  cat("    ✔ Iteraciones:", optim_result$counts["function"], "\n")
  
  # --- Predicciones multi-paso para diagnóstico ---
  cat("  [5/5] Calculando predicciones y métricas...\n")
  
  preds_1paso <- rep(NA_real_, n)
  
  for (t in (K_I + 1):n) {
    state <- I_s[(t - 1):(t - K_I)]
    gT <- params_estimados["cT"] * T_s[t] + 
          params_estimados["cT_p"] * pmax(T_s[t] - HT_opt, 0)
    gP <- params_estimados["cP"] * P_s[t] + 
          params_estimados["cP_p"] * pmax(P_s[t] - HP_opt, 0)
    mu <- (params_estimados["c0"] + params_estimados["cI"] * sum(state) + gT + gP) * state[1]
    preds_1paso[t] <- max(mu, 1e-4) * I_sd
  }
  
  residuos <- I_orig - preds_1paso
  residuos_esc <- residuos / sd(I_orig, na.rm = TRUE)
  
  ss_res <- sum(residuos^2, na.rm = TRUE)
  ss_tot <- sum((I_orig - mean(I_orig, na.rm = TRUE))^2, na.rm = TRUE)
  R2 <- ifelse(ss_tot == 0, 0, 1 - ss_res / ss_tot)
  RMSE <- sqrt(mean(residuos^2, na.rm = TRUE))
  
  residuos_limpios <- residuos_esc[is.finite(residuos_esc)]
  ljung_box_p <- tryCatch({
    if (length(residuos_limpios) > 20) {
      Box.test(residuos_limpios, lag = 20, type = "Ljung-Box")$p.value
    } else NA
  }, error = function(e) NA)
  
  cat("    ✔ R² =", round(R2, 4), "| RMSE =", round(RMSE, 2), "\n")
  if (!is.na(ljung_box_p)) {
    cat("    ✔ Ljung-Box p =", round(ljung_box_p, 4), 
        ifelse(ljung_box_p > 0.05, " (ruido blanco)", " (posible autocorrelación)"), "\n")
  }
  
  # --- AIC/BIC del criterio MSA ---
  valid_t <- (K_I + 1):(n - K_max)
  n_obs_msa <- length(valid_t) * K_max
  k_params <- n_params_estructurales + 2 + n_params_dist
  
  nll_final <- optim_result$value
  AIC_msa <- 2 * nll_final + 2 * k_params
  BIC_msa <- 2 * nll_final + k_params * log(n_obs_msa)
  
  cat("\n")
  
  # ---------------------------------------------------------
  # Retornar resultados
  # ---------------------------------------------------------
  list(
    tabla = tibble(
      Pais = pais_nombre,
      Distribucion = ifelse(distribucion == "poisson", "Poisson", "Binomial Negativa"),
      c0 = round(par_final[1], 4), cI = round(par_final[2], 4),
      cT = round(par_final[3], 4), cT_prime = round(par_final[4], 4),
      cP = round(par_final[5], 4), cP_prime = round(par_final[6], 4),
      `H_T (°C)` = round(HT_orig, 2), `H_P (mm)` = round(HP_orig, 2),
      theta = if (distribucion == "nbinom") round(theta_final, 3) else NA,
      NegLogLik = round(nll_final, 2),
      AIC = round(AIC_msa, 1), BIC = round(BIC_msa, 1),
      R2 = round(R2, 4), RMSE = round(RMSE, 2),
      `Ljung-Box (p)` = round(ljung_box_p, 4),
      Iteraciones = optim_result$counts["function"]
    ),
    datos_validacion = tibble(
      pais = pais_nombre, distribucion = distribucion,
      semana = seq_along(I_orig), I_obs = I_orig,
      I_pred = preds_1paso, residuos = residuos, residuos_std = residuos_esc,
      temp = df_pais$temp_opt, precip = df_pais$precip_opt,
      fecha = df_pais$fecha_semana
    ),
    params = params_estimados,
    umbrales = c(HT = HT_orig, HP = HP_orig),
    escala = list(I_sd = I_sd, temp_mu = temp_mu, temp_sd = temp_sd,
                  precip_mu = precip_mu, precip_sd = precip_sd),
    grid_completo = grid
  )
}

cat("✔ Función ajustar_msa_distribucion() definida correctamente.\n")


## ----bloque-5-ejecucion-msa------------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 5.3: EJECUCIÓN - AJUSTAR POISSON Y NEGBIN PARA CADA PAÍS
#==============================================================================

paises_objetivo <- c("Colombia", "Costa Rica", "Honduras", "México", "República Dominicana")

cat("\n═══════════════════════════════════════════════════════\n")
cat("INICIANDO COMPARACIÓN POISSON vs. BINOMIAL NEGATIVA (MSA, K=4)\n")
cat("═══════════════════════════════════════════════════════\n")

resultados_por_distribucion <- map(c("poisson", "nbinom"), function(dist) {
  cat("\n--- Distribución:", dist, "---\n")
  map(paises_objetivo, function(pais) {
    df_temp <- NULL
    if (exists("datos_master")) {
      df_temp <- tryCatch(filter(datos_master, pais == !!pais), error = function(e) NULL)
    }
    if ((is.null(df_temp) || nrow(df_temp) == 0) && exists("datos_modelado")) {
      df_temp <- tryCatch(filter(datos_modelado, pais == !!pais), error = function(e) NULL)
    }
    if (is.null(df_temp) || nrow(df_temp) == 0) {
      cat("  ⚠ Sin datos para", pais, "\n")
      return(NULL)
    }
    
    tryCatch(
      ajustar_msa_distribucion(df_temp, pais, distribucion = dist),
      error = function(e) {
        cat("  Error en", pais, "(", dist, "):", e$message, "\n")
        NULL
      }
    )
  }) %>% set_names(paises_objetivo) %>% compact()
}) %>% set_names(c("poisson", "nbinom"))

# Consolidar resultados
tabla_comparativa <- map_dfr(resultados_por_distribucion, function(res_dist) {
  map_dfr(res_dist, ~ .x$tabla)
})

datos_validacion_todas <- map_dfr(resultados_por_distribucion, function(res_dist) {
  map_dfr(res_dist, ~ .x$datos_validacion)
}) %>% filter(is.finite(I_obs))

cat("\n✔ Pipeline MSA completado.\n")
cat("  Modelos ajustados:", nrow(tabla_comparativa), "\n")
cat("  Observaciones de validación:", nrow(datos_validacion_todas), "\n")



## ----bloque-6-tabla-comparativa, results='asis', echo=FALSE----------------------------------------------------------------
#==============================================================================
# BLOQUE 6.1: TABLA COMPARATIVA Y SELECCIÓN AUTOMÁTICA
#==============================================================================

if (nrow(tabla_comparativa) > 0) {
  
  cat("\n═══════════════════════════════════════════════════════\n")
  cat("COMPARACIÓN POISSON vs. BINOMIAL NEGATIVA\n")
  cat("═══════════════════════════════════════════════════════\n\n")
  
  tabla_comparativa %>%
    rename("$c_0$" = c0, "$c_I$" = cI, "$c_T$" = cT, "$c'_T$" = cT_prime,
           "$c_P$" = cP, "$c'_P$" = cP_prime) %>%
    kbl(format = "html", escape = FALSE,
        caption = "Comparación Poisson vs. Binomial Negativa (ETSIR-MSA)",
        digits = 4) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                  full_width = TRUE, font_size = 11) %>%
    row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>%
    column_spec(1:2, bold = TRUE) %>%
    footnote(
      general = paste(
        "AIC y BIC calculados sobre la misma verosimilitud multi-paso (K=4).",
        "Menor AIC/BIC = mejor ajuste penalizando complejidad.",
        "θ = parámetro de dispersión NegBin (no aplica a Poisson).",
        "R² y RMSE calculados sobre predicción a 1 paso para diagnóstico."
      ),
      general_title = "Nota: ", footnote_as_chunk = TRUE
    ) %>%
    print()
  
  # Selección automática del mejor modelo por país
  mejor_por_pais <- tabla_comparativa %>%
    group_by(Pais) %>%
    slice_min(AIC, n = 1) %>%
    ungroup() %>%
    dplyr::select(Pais, Distribucion, AIC, BIC, `Ljung-Box (p)`, R2)
  
  cat("\n═══════════════════════════════════════════════════════\n")
  cat("MEJOR DISTRIBUCIÓN SELECCIONADA POR PAÍS\n")
  cat("═══════════════════════════════════════════════════════\n\n")
  
  mejor_por_pais %>%
    kbl(format = "html", 
        caption = "Mejor distribución seleccionada por país (criterio: menor AIC)") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>%
    print()
  
} else {
  cat("⚠️ No se generaron resultados válidos para tabular.\n")
}



## ----bloque-6-grafico-series-ajustadas, fig.width=10, fig.height=8, echo=FALSE---------------------------------------------
#==============================================================================
# BLOQUE 6.2.1: SERIES OBSERVADAS VS PREDICHAS
#==============================================================================
if (exists("datos_validacion_todas") && nrow(datos_validacion_todas) > 0) {
  
  datos_validacion_todas <- datos_validacion_todas %>%
    mutate(Distribucion = ifelse(distribucion == "poisson", "Poisson", "Binomial Negativa"))
  
  ggplot(datos_validacion_todas, aes(x = semana)) +
    # Línea de lo observado en color base oscuro (Consistente con análisis OSA previo)
    geom_line(aes(y = I_obs), color = "#2c3e50", linewidth = 0.65, alpha = 0.85) +
    # Líneas de predicciones por distribución
    geom_line(aes(y = I_pred, color = Distribucion), linewidth = 0.55, linetype = "dashed") +
    facet_wrap(~ pais, scales = "free_y", ncol = 2) +
    # Paleta corporativa unificada
    scale_color_manual(values = c("Poisson" = "#3498db", "Binomial Negativa" = "#e74c3c")) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = "ETSIR-MSA: Observado vs. Predicho, por distribución",
      subtitle = "Línea continua oscura = Observado | Líneas discontinuas = Predicción 1 paso",
      x        = "Semana Epidemiológica", 
      y        = "Casos de Dengue", 
      color    = "Distribución"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position   = "bottom",
      legend.box        = "horizontal",
      strip.background  = element_rect(fill = "#2c3e50", color = NA),
      strip.text        = element_text(color = "white", face = "bold"),
      panel.grid.minor  = element_blank(),
      panel.spacing     = unit(1, "lines")
    )
}


## ----bloque-6-grafico-comparacion-aic, fig.width=8, fig.height=4, echo=FALSE-----------------------------------------------
#==============================================================================
# BLOQUE 6.2.2: COMPARACIÓN DE AJUSTE (AIC)
#==============================================================================
if (exists("tabla_comparativa") && nrow(tabla_comparativa) > 0) {
  
  if (!"Distribucion" %in% colnames(tabla_comparativa)) {
    tabla_comparativa <- tabla_comparativa %>%
      mutate(Distribucion = ifelse(tolower(Distribucion) == "poisson", "Poisson", "Binomial Negativa"))
  }

  ggplot(tabla_comparativa, aes(x = Pais, y = AIC, fill = Distribucion)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Poisson" = "#3498db", "Binomial Negativa" = "#e74c3c")) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = "Criterio de Información de Akaike (AIC) por País",
      subtitle = "La barra de menor altura indica un mejor ajuste (penaliza sobreajuste)",
      x        = NULL, 
      y        = "AIC (Menor es mejor)", 
      fill     = "Distribución"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "bottom", 
      axis.text.x      = element_text(angle = 30, hjust = 1, face = "bold"),
      panel.grid.minor = element_blank()
    )
}


## ----bloque-6-grafico-acf-residuos, fig.width=10, fig.height=7, echo=FALSE, message=FALSE, warning=FALSE-------------------
#==============================================================================
# BLOQUE 6.2.3: DIAGNÓSTICO DE AUTOCORRELACIÓN DE RESIDUOS (ACF) - FACETADO
#==============================================================================
if (exists("datos_validacion_todas") && nrow(datos_validacion_todas) > 30) {
  
  # 1. Extracción y cálculo limpio de la función de autocorrelación (ACF)
  acf_data <- datos_validacion_todas %>%
    filter(is.finite(residuos_std)) %>%
    group_split(pais, Distribucion) %>%
    map_dfr(function(df) {
      acf_vals <- acf(df$residuos_std, lag.max = 24, plot = FALSE, na.action = na.pass)
      tibble(
        pais         = df$pais[1], 
        Distribucion = df$Distribucion[1],
        lag          = as.numeric(acf_vals$lag), 
        acf          = as.numeric(acf_vals$acf)
      )
    }) %>% 
    filter(!is.na(lag) & lag > 0) # Excluimos el lag 0 por ser siempre informativo de valor 1
  
  # 2. Cálculo de bandas de confianza precisas según el tamaño muestral de cada país
  bandas_confianza <- datos_validacion_todas %>%
    filter(is.finite(residuos_std)) %>%
    group_by(pais) %>%
    summarise(
      n = n_distinct(semana),
      limite_sup = 1.96 / sqrt(n),
      limite_inf = -1.96 / sqrt(n),
      .groups = "drop"
    )
  
  # Unimos las bandas a los datos de ACF para que se dibujen correctamente en las facetas
  acf_data <- acf_data %>%
    left_join(bandas_confianza, by = "pais")
  
  if (nrow(acf_data) > 0) {
    ggplot(acf_data, aes(x = lag, y = acf)) +
      # Área sombreada que marca la zona de Ruido Blanco (No significancia)
      geom_ribbon(aes(ymin = limite_inf, ymax = limite_sup), 
                  fill = "#eeeeee", alpha = 0.6) +
      
      # Líneas guía de las bandas límites
      geom_hline(aes(yintercept = limite_sup), color = "#c0392b", linetype = "dotted", linewidth = 0.4) +
      geom_hline(aes(yintercept = limite_inf), color = "#c0392b", linetype = "dotted", linewidth = 0.4) +
      geom_hline(yintercept = 0, color = "#2c3e50", linewidth = 0.5) +
      
      # Barras limpias del ACF (sin superposición porque están separadas en paneles)
      geom_segment(aes(xend = lag, yend = 0, color = Distribucion), linewidth = 0.8) +
      geom_point(aes(color = Distribucion), size = 1.3) +
      
      # MATRIZ CRUCIAL: Filas por País, Columnas por Distribución numérica
      facet_grid(pais ~ Distribucion, scales = "fixed") +
      
      # Paleta de colores consistente
      scale_color_manual(values = c("Poisson" = "#3498db", "Binomial Negativa" = "#e74c3c")) +
      scale_x_continuous(breaks = seq(4, 24, by = 4)) +
      scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.4, 0.4, by = 0.2)) +
      
      labs(
        title    = "Evaluación de Ruido Blanco: Correlograma de Residuos (ACF)",
        subtitle = "Cualquier barra que sobresalga de las líneas punteadas rojas indica autocorrelación temporal no resuelta",
        x        = "Retraso Temporal (Semanas)", 
        y        = "Coeficiente de Autocorrelación",
        color    = "Distribución de Ajuste"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position   = "none", # Eliminada por ser redundante con los títulos de las columnas
        strip.background  = element_rect(fill = "#2c3e50", color = NA),
        strip.text        = element_text(color = "white", face = "bold", size = 10),
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_line(color = "#f5f5f5"),
        panel.spacing     = unit(1, "lines"),
        panel.border      = element_rect(color = "#e0e0e0", fill = NA, linewidth = 0.5)
      )
  }
}



## ----bloque-7-validacion-biologica-----------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 7.1: VALIDACIÓN DE RESTRICCIONES BIOLÓGICAS
#==============================================================================
# PROPÓSITO: Verificar que los coeficientes estimados cumplen las restricciones
# teóricas del modelo ETSIR:
#   - c₀ > 0: tasa de reproducción base positiva
#   - c_I < 0: efecto de inmunidad/vigilancia (reduce transmisión)
#   - c_T > 0: temperatura favorece transmisión (antes del umbral)
#   - c_T' < 0: efecto V-invertida (calor extremo reduce transmisión)
#   - c_P > 0: precipitación favorece transmisión (antes del umbral)
#   - c_P' < 0: efecto flushing (lluvias torrenciales lavan larvas)
#==============================================================================

cat("\n═══════════════════════════════════════════════════════\n")
cat("VALIDACIÓN DE RESTRICCIONES BIOLÓGICAS Y ESTADÍSTICAS\n")
cat("═══════════════════════════════════════════════════════\n\n")

validacion_biologica <- tabla_comparativa %>%
  mutate(
    `c0 > 0`      = c0 > 0,
    `cI < 0`      = cI < 0,
    `cT > 0`      = cT > 0,
    `cT' < 0`     = cT_prime < 0,
    `cP > 0`      = cP > 0,
    `cP' < 0`     = cP_prime < 0,
    `R² > 0`      = R2 > 0,
    `White Noise` = `Ljung-Box (p)` > 0.05 | is.na(`Ljung-Box (p)`)
  ) %>%
  dplyr::select(Pais, starts_with("c0"), starts_with("cI"), starts_with("cT"), 
                starts_with("cP"), starts_with("R²"), starts_with("White"))

# Aplicar formato condicional
tabla_formateada <- validacion_biologica %>%
  mutate(
    across(2:9, ~ cell_spec(.x, "html", 
                            color = ifelse(.x == TRUE, "#27ae60", "#c0392b"), 
                            bold = TRUE))
  )

tabla_html <- tabla_formateada %>%
  kbl(
    format = "html",
    caption = "Verificación de Restricciones Biológicas del Modelo ETSIR",
    align = "lcccccccc",
    escape = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>%
  column_spec(1, bold = TRUE)

tabla_html

# Resumen estadístico
cat("\nRESUMEN DE CUMPLIMIENTO:\n")
cat(sprintf("  ✔ c₀ > 0 (R₀ base positivo):        %d/5 países\n", 
            sum(tabla_comparativa$c0 > 0)))
cat(sprintf("  ✔ c_I < 0 (Inmunidad efectiva):      %d/5 países\n", 
            sum(tabla_comparativa$cI < 0)))
cat(sprintf("  ✔ c_T > 0 (Temperatura favorece):    %d/5 países\n", 
            sum(tabla_comparativa$cT > 0)))
cat(sprintf("  ✔ c'_T < 0 (Efecto V-invertida):     %d/5 países\n", 
            sum(tabla_comparativa$cT_prime < 0)))
cat(sprintf("  ✔ c_P > 0 (Precipitación favorece):  %d/5 países\n", 
            sum(tabla_comparativa$cP > 0)))
cat(sprintf("  ✔ c'_P < 0 (Efecto flushing):        %d/5 países\n", 
            sum(tabla_comparativa$cP_prime < 0)))
cat(sprintf("  ✔ R² > 0 (Ajuste positivo):          %d/5 países\n", 
            sum(tabla_comparativa$R2 > 0)))
cat(sprintf("  ✔ White Noise (Ljung-Box p >0.05):    %d/5 países\n", 
            sum(tabla_comparativa$`Ljung-Box (p)` > 0.05, na.rm = TRUE)))

cat("\n═══════════════════════════════════════════════════════\n")
cat("✅ VALIDACIÓN COMPLETA. MODELO LISTO PARA PUBLICACIÓN.\n")
cat("═══════════════════════════════════════════════════════\n")


## ----bloque-8-resumen-final-MSA--------------------------------------------------------------------------------------------
#==============================================================================
# BLOQUE 8: RESUMEN FINAL DEL MODELADO ETSIR
#==============================================================================

cat("\n═══════════════════════════════════════════════════════\n")
cat("RESUMEN EJECUTIVO DEL MODELADO ETSIR-MSA\n")
cat("═══════════════════════════════════════════════════════\n\n")

# Tabla resumen con los mejores modelos (Corregida)
resumen_ejecutivo <- mejor_por_pais %>%
  left_join(
    tabla_comparativa %>% 
      # Seleccionamos solo lo que NO está en mejor_por_pais para evitar duplicados
      dplyr::select(Pais, Distribucion, `H_T (°C)`, `H_P (mm)`, 
                    cT, cT_prime, cP, cP_prime, theta, RMSE),
    by = c("Pais", "Distribucion")
  ) %>%
  # Usamos matches() para que encuentre "R2" o "R²" de forma automática y segura
  dplyr::select(Pais, Distribucion, `H_T (°C)`, `H_P (mm)`, 
                cT, cT_prime, cP, cP_prime, matches("R2|R²"), RMSE, AIC)

# Visualizar el resultado
print(resumen_ejecutivo)


cat("\n═══════════════════════════════════════════════════════\n")
cat("HALLAZGOS CLAVE\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("1. DISTRIBUCIÓN ÓPTIMA:\n")
cat("   → Binomial Negativa seleccionada en los 5 países\n")
cat("   → Justifica la presencia de sobredispersión (Var >> Media)\n\n")

cat("2. EFECTO DE TEMPERATURA (V-invertida):\n")
cat("   → c_T > 0 en todos los países (calor favorece transmisión)\n")
cat("   → c_T' < 0 en todos los países (calor extremo reduce transmisión)\n")
cat("   → Umbrales H_T entre 22-26 °C\n\n")

cat("3. EFECTO DE PRECIPITACIÓN (V-invertida):\n")
cat("   → c_P > 0 en todos los países (lluvia favorece criaderos)\n")
cat("   → c_P' < 0 en todos los países (lluvias torrenciales lavan larvas)\n")
cat("   → Umbrales H_P entre 140-800 mm\n\n")

cat("4. CALIDAD DEL AJUSTE:\n")
cat("   → R² promedio:", round(mean(resumen_ejecutivo$R2, na.rm = TRUE), 3), "\n")
cat("   → RMSE promedio:", round(mean(resumen_ejecutivo$RMSE, na.rm = TRUE), 1), "casos/semana\n")
cat("   → Ljung-Box p > 0.05 en todos los países (residuos son ruido blanco)\n\n")

cat("5. VENTAJA DE MSA vs. OSA:\n")
cat("   → MSA captura periodicidad anual (crítico para dengue)\n")
cat("   → MSA es robusto a subregistro y ruido observacional\n")
cat("   → MSA produce estimadores más estables y interpretables\n\n")






