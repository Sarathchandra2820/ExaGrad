#!/bin/bash
# ============================================================
#  ExaGrad — Compile & Run Direct RHF SCF
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- defaults ------------------------------------------------
XYZ="geometry/pyridine.xyz"
BASIS="sto-6g"
METHOD="direct"

# --- parse optional arguments --------------------------------
usage() {
    echo "Usage: $0 [-g geometry.xyz] [-b basis_set] [-m direct|df] [-c] [-h]"
    echo ""
    echo "  -g FILE    Path to .xyz geometry file  (default: $XYZ)"
    echo "  -b BASIS   Basis set name              (default: $BASIS)"
    echo "  -m METHOD  Fock builder: direct|df     (default: $METHOD)"
    echo "  -c         Clean build (make clean first)"
    echo "  -h         Show this help message"
    exit 0
}

CLEAN=false
while getopts "g:b:m:ch" opt; do
    case $opt in
        g) XYZ="$OPTARG" ;;
        b) BASIS="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        c) CLEAN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

METHOD_LC="$(printf '%s' "$METHOD" | tr '[:upper:]' '[:lower:]')"
if [[ "$METHOD_LC" != "direct" && "$METHOD_LC" != "df" ]]; then
    echo "ERROR: Invalid -m value '$METHOD'. Use 'direct' or 'df'."
    exit 2
fi

# --- clean (optional) ----------------------------------------
if $CLEAN; then
    echo "==> Cleaning previous build …"
    make clean
    echo ""
fi

# --- generate integrals with PySCF / libcint -----------------
echo "==> Generating integrals for $XYZ with basis $BASIS …"
python3 python/export_cint_env.py "$XYZ" "$BASIS"
echo ""

# --- compile -------------------------------------------------
echo "==> Compiling Fortran sources …"
make -j rhf_main
echo ""

# --- run -----------------------------------------------------
echo "==> Running RHF SCF …"
echo "------------------------------------------------------------"
LOG_FILE="$(mktemp -t exagrad_scf.XXXXXX)"
EXAGRAD_FOCK_BUILDER="$METHOD_LC" ./rhf_main | tee "$LOG_FILE"
echo "------------------------------------------------------------"

echo "==> Running PySCF reference and comparison …"
set +e
python3 python/compare_scf_energy.py \
    --xyz "$XYZ" \
    --basis "$BASIS" \
    --fortran-log "$LOG_FILE"
CMP_RC=$?
set -e

if [ "$CMP_RC" -eq 1 ]; then
    echo "==> Comparison warning: energies differ beyond tolerance."
elif [ "$CMP_RC" -ne 0 ]; then
    echo "==> Comparison failed (code $CMP_RC)."
fi

rm -f "$LOG_FILE"
echo "==> Done."
