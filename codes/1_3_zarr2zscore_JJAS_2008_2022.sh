#!/usr/bin/env bash
set -euo pipefail

# ===========================
# CONFIGURACIÓN
# ===========================

# Variables dentro del Zarr que quieres procesar
VARS=(lai smi lst_day lst_night ndvi rh t2m wind_speed tp)

STORE_PATH="/diskonfire/MESOGEOS"
OUT_DIR="/diskonfire/MESOGEOS"

# Ruta al venv donde están xarray/dask
VENV_ACTIVATE="/diskonfire/MESOGEOS/.venv/bin/activate"

# Directorio donde está este script (para llamar al .py con ruta absoluta)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ===========================
# INICIO
# ===========================

echo "=== Activando entorno Python: ${VENV_ACTIVATE} ==="
source "${VENV_ACTIVATE}"

for VAR in "${VARS[@]}"; do
  echo
  echo "==============================="
  echo "Procesando variable: ${VAR}"
  echo "==============================="

  IN_NC="${OUT_DIR}/${VAR}_JJAS_2008_2022.nc"
  CLIM_NC="${OUT_DIR}/${VAR}_JJAS_clim_2008_2022.nc"
  STD_NC="${OUT_DIR}/${VAR}_JJAS_std_2008_2022.nc"
  Z_NC="${OUT_DIR}/${VAR}_JJAS_zscore_2008_2022.nc"

  # ---------------------------
  # 1) Zarr -> NetCDF (JJAS 2008–2022)
  # ---------------------------
  if [[ -f "${IN_NC}" ]]; then
    echo "   [INFO] Ya existe ${IN_NC}, salto export desde Zarr."
  else
    echo "   [STEP] Exportando ${VAR} JJAS 2008–2022 a NetCDF..."
    python "${SCRIPT_DIR}/export_var_JJAS_2008_2022_to_nc.py" \
      --var "${VAR}" \
      --store "${STORE_PATH}" \
      --outdir "${OUT_DIR}"
    echo "   [OK] Creado ${IN_NC}"
  fi

  # Por seguridad, borramos salidas parciales si las hubiera
  rm -f "${CLIM_NC}" "${STD_NC}" "${Z_NC}" || true

  # ---------------------------
  # 2) CDO: climatología en JJAS (media temporal)
  #    -> timmean = media sobre todos los días JJAS 2008–2022
  # ---------------------------
  echo "   [STEP] Calculando timmean (climatología JJAS) con CDO (low memory)..."
  cdo -L -b f32 timmean "${IN_NC}" "${CLIM_NC}"

  # ---------------------------
  # 3) CDO: std temporal en JJAS
  # ---------------------------
  echo "   [STEP] Calculando timstd (std JJAS) con CDO (low memory)..."
  cdo -L -b f32 timstd "${IN_NC}" "${STD_NC}"

  # ---------------------------
  # 4) CDO: z-score = (X - clim) / std
  # ---------------------------
  echo "   [STEP] Calculando z-score diario (div -sub, low memory)..."
  cdo -L -b f32 div \
    -sub "${IN_NC}" "${CLIM_NC}" \
    "${STD_NC}" \
    "${Z_NC}"

  echo "   [OK] z-score guardado en ${Z_NC}"
done

echo
echo "=== Todo completado para: ${VARS[*]} ==="
