#!/bin/bash
# ============================================================
#  ExaGrad — Compile & Run ROSE-RHF (supersystem + fragments)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- defaults ------------------------------------------------
XYZ=""
BASIS="cc-pVDZ"
METHOD="true_df"
DATA_DIR="geometry/pyridine_rose"
USER_SET_GEOMETRY=false
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
    echo "Usage: $0 [-g geometry.xyz] [-b basis_set] [-m direct|cholesky|block_cholesky|true_df] [-a auxbasis] [-d data_dir] [-t nthreads] [-c] [-h]"
    echo ""
    echo "  -g FILE    Path to supersystem .xyz file     (default: DATA_DIR/supra_mol.xyz if present)"
    echo "  -b BASIS   Basis set name                    (default: $BASIS)"
    echo "  -m METHOD  Fock builder: direct|cholesky|block_cholesky|true_df  (default: $METHOD)"
    echo "  -a AUX     Auxiliary basis for true_df       (default: $AUXBASIS)"
    echo "  -d DIR     Molecule data root directory      (default: $DATA_DIR)"
    echo "  -t N       Number of OpenMP/BLAS threads     (default: $NTHREADS)"
    echo "  -c         Clean build (make clean first)"
    echo "  -h         Show this help message"
    exit 0
}

resolve_geometry() {
    local data_dir="$1"

    if $USER_SET_GEOMETRY; then
        printf '%s\n' "$XYZ"
        return 0
    fi

    if [[ -f "$data_dir/supra_mol.xyz" ]]; then
        printf '%s\n' "$data_dir/supra_mol.xyz"
        return 0
    fi

    local candidates=()
    local xyz_path
    shopt -s nullglob
    for xyz_path in "$data_dir"/*.xyz; do
        [[ "$(basename "$xyz_path")" == frag_* ]] && continue
        candidates+=("$xyz_path")
    done
    shopt -u nullglob

    if [[ ${#candidates[@]} -eq 1 ]]; then
        printf '%s\n' "${candidates[0]}"
        return 0
    fi

    if [[ ${#candidates[@]} -gt 1 ]]; then
        echo "ERROR: Multiple supersystem XYZ files found in '$data_dir'." >&2
        echo "Set one as supra_mol.xyz or pass -g explicitly." >&2
        exit 2
    fi

    echo "ERROR: No supersystem XYZ file found in '$data_dir'." >&2
    echo "Expected supra_mol.xyz or pass -g explicitly." >&2
    exit 2
}

CLEAN=false
while getopts "g:b:m:a:d:t:ch" opt; do
    case $opt in
        g) XYZ="$OPTARG"; USER_SET_GEOMETRY=true ;;
        b) BASIS="$OPTARG" ;;
        m) METHOD="$OPTARG" ;;
        a) AUXBASIS="$OPTARG" ;;
        d) DATA_DIR="$OPTARG" ;;
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

XYZ="$(resolve_geometry "$DATA_DIR")"

if ! $USER_SET_GEOMETRY && [[ "$(basename "$XYZ")" == "supra_mol.xyz" ]]; then
    echo "==> Detected supramolecule cue: using $XYZ"
fi

if [[ ! -f "$XYZ" ]]; then
    echo "ERROR: Geometry file '$XYZ' was not found."
    exit 2
fi

# --- clean (optional) ----------------------------------------
if $CLEAN; then
    echo "==> Cleaning previous build …"
    make clean
    echo ""
fi

# --- generate integrals for supersystem + fragments ----------
echo "==> Generating integrals for $XYZ with basis $BASIS …"
echo "    (supersystem + any frag_*.xyz fragments in $DATA_DIR)"
mkdir -p "$DATA_DIR/ints"
python3 python/export_cint_env.py "$XYZ" "$BASIS" "$DATA_DIR"

if [[ "$METHOD_LC" == "true_df" ]]; then
    echo "==> Exporting true DF metadata (supersystem + fragments) …"
    python3 python/export_df_metadata.py "$XYZ" "$BASIS" "$AUXBASIS" "$DATA_DIR"
fi
echo ""

# --- compile -------------------------------------------------
echo "==> Compiling Fortran sources …"
make rhf_rose_main
echo ""

# --- run ROSE-RHF (single process, loads everything) ---------
echo "==> Running ROSE-RHF …"
echo "  Fock builder : $METHOD_LC"
echo "  Threads      : $NTHREADS"
echo "  Data dir     : $DATA_DIR"
echo "------------------------------------------------------------"
EXAGRAD_FOCK_BUILDER="$METHOD_LC" \
EXAGRAD_MOL_DIR="$DATA_DIR" \
OMP_NUM_THREADS="$NTHREADS" \
OPENBLAS_NUM_THREADS="$NTHREADS" \
./rhf_rose_main
echo "------------------------------------------------------------"

echo "==> Done."
