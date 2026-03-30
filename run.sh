#!/bin/bash
# ============================================================
#  ExaGrad — Compile & Run Direct RHF SCF
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- defaults ------------------------------------------------
XYZ="geometry/pyridine_dimer.xyz"
BASIS="cc-pVDZ"
METHOD="true_df"
if command -v nproc >/dev/null 2>&1; then
    NTHREADS="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    NTHREADS="$(sysctl -n hw.ncpu)"
else
    NTHREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi
AUXBASIS="weigend"

# --- parse optional arguments --------------------------------
usage() {
    echo "Usage: $0 [-g geometry.xyz] [-b basis_set] [-m direct|cholesky|block_cholesky|true_df] [-a auxbasis] [-t nthreads] [-c] [-h]"
    echo ""
    echo "  -g FILE    Path to .xyz geometry file        (default: $XYZ)"
    echo "  -b BASIS   Basis set name                    (default: $BASIS)"
    echo "  -m METHOD  Fock builder: direct|cholesky|block_cholesky|true_df  (default: $METHOD)"
    echo "  -a AUX     Auxiliary basis for true_df         (default: $AUXBASIS)"
    echo "  -t N       Number of OpenMP/BLAS threads     (default: $NTHREADS)"
    echo "  -c         Clean build (make clean first)"
    echo "  -h         Show this help message"
    exit 0
}

CLEAN=false
while getopts "g:b:m:a:t:ch" opt; do
    case $opt in
        g) XYZ="$OPTARG" ;;
        b) BASIS="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        a) AUXBASIS="$OPTARG" ;;
        t) NTHREADS="$OPTARG" ;;
        c) CLEAN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if ! [[ "$NTHREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: -t value must be a positive integer (got '$NTHREADS')."
    exit 2
fi

METHOD_LC="$(printf '%s' "$METHOD" | tr '[:upper:]' '[:lower:]')"
if [[ "$METHOD_LC" == "df" ]]; then
    METHOD_LC="cholesky"
elif [[ "$METHOD_LC" == "block_df" || "$METHOD_LC" == "df_block" || "$METHOD_LC" == "blocked_df" ]]; then
    METHOD_LC="block_cholesky"
fi

if [[ "$METHOD_LC" != "direct" && "$METHOD_LC" != "cholesky" && "$METHOD_LC" != "block_cholesky" && "$METHOD_LC" != "true_df" ]]; then
    echo "ERROR: Invalid -m value '$METHOD'. Use 'direct', 'cholesky', 'block_cholesky', or 'true_df'."
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

if [[ "$METHOD_LC" == "true_df" ]]; then
    echo "==> Exporting true DF metadata (combined env only) …"
    python3 python/export_df_metadata.py "$XYZ" "$BASIS" "$AUXBASIS"
fi
echo ""

# --- compile -------------------------------------------------
echo "==> Compiling Fortran sources …"
make -j rhf_main
echo ""

# --- run -----------------------------------------------------
echo "==> Running RHF SCF …"
echo "  Fock builder : $METHOD_LC"
echo "  Threads      : $NTHREADS"
echo "------------------------------------------------------------"
LOG_FILE="$(mktemp -t exagrad_scf.XXXXXX)"
EXAGRAD_FOCK_BUILDER="$METHOD_LC" OMP_NUM_THREADS="$NTHREADS" OPENBLAS_NUM_THREADS="$NTHREADS" ./rhf_main | tee "$LOG_FILE"
echo "------------------------------------------------------------"

echo "==> Running PySCF reference and comparison …"
set +e
python3 python/compare_scf_energy.py \
    --xyz "$XYZ" \
    --basis "$BASIS" \
    --fortran-method "$METHOD_LC" \
    --auxbasis "$AUXBASIS" \
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
