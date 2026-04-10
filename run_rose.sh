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
    DEFAULT_THREADS="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    DEFAULT_THREADS="$(sysctl -n hw.ncpu)"
else
    DEFAULT_THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
fi
OMP_THREADS="$DEFAULT_THREADS"
OPENBLAS_THREADS="$DEFAULT_THREADS"
AUXBASIS="weigend"
USER_SET_OMP_THREADS=false
USER_SET_OPENBLAS_THREADS=false

# --- parse optional arguments --------------------------------
usage() {
    echo "Usage: $0 [-g geometry.xyz] [-b basis_set] [-m direct|cholesky|block_cholesky|true_df|blocked_true_df] [-a auxbasis] [-d data_dir] [-t nthreads] [-t1 openblas_threads] [-t2 openmp_threads] [-o fragment|fragment_recanon|energy] [-c] [-h]"
    echo ""
    echo "  -g FILE    Path to supersystem .xyz file     (default: DATA_DIR/supra_mol.xyz if present)"
    echo "  -b BASIS   Basis set name                    (default: $BASIS)"
    echo "  -m METHOD  Fock builder: direct|cholesky|block_cholesky|true_df|blocked_true_df  (default: $METHOD)"
    echo "  -a AUX     Auxiliary basis for true_df       (default: $AUXBASIS)"
    echo "  -d DIR     Molecule data root directory      (default: $DATA_DIR)"
    echo "  -t N       Set both OpenBLAS and OpenMP threads to N (default: $DEFAULT_THREADS)"
    echo "  -t1 N      Set OpenBLAS threads only         (default: $DEFAULT_THREADS)"
    echo "  -t2 N      Set OpenMP threads only           (default: $DEFAULT_THREADS)"
    echo "  -o MODE    LMO ordering: fragment|fragment_recanon|energy (default: fragment)"
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
ROSE_ORDER="fragment"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g)
            XYZ="$2"; USER_SET_GEOMETRY=true; shift 2 ;;
        -b)
            BASIS="$2"; shift 2 ;;
        -m)
            METHOD="$2"; shift 2 ;;
        -a)
            AUXBASIS="$2"; shift 2 ;;
        -d)
            DATA_DIR="$2"; shift 2 ;;
        -t)
            OMP_THREADS="$2"
            OPENBLAS_THREADS="$2"
            USER_SET_OMP_THREADS=true
            USER_SET_OPENBLAS_THREADS=true
            shift 2 ;;
        -t1)
            OPENBLAS_THREADS="$2"
            USER_SET_OPENBLAS_THREADS=true
            shift 2 ;;
        -t2)
            OMP_THREADS="$2"
            USER_SET_OMP_THREADS=true
            shift 2 ;;
        -o)
            ROSE_ORDER="$2"; shift 2 ;;
        -c)
            CLEAN=true; shift ;;
        -h)
            usage ;;
        --)
            shift; break ;;
        *)
            echo "ERROR: Unknown option '$1'"
            usage ;;
    esac
done

if ! [[ "$OPENBLAS_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: -t/-t1 value must be a positive integer (got '$OPENBLAS_THREADS')."
    exit 2
fi

if ! [[ "$OMP_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: -t/-t2 value must be a positive integer (got '$OMP_THREADS')."
    exit 2
fi

METHOD_LC="$(printf '%s' "$METHOD" | tr '[:upper:]' '[:lower:]')"
if [[ "$METHOD_LC" == "df" ]]; then
    METHOD_LC="cholesky"
elif [[ "$METHOD_LC" == "block_df" || "$METHOD_LC" == "df_block" || "$METHOD_LC" == "blocked_df" ]]; then
    METHOD_LC="block_cholesky"
fi

if [[ "$METHOD_LC" != "direct" && "$METHOD_LC" != "cholesky" && "$METHOD_LC" != "block_cholesky" && "$METHOD_LC" != "true_df" && "$METHOD_LC" != "blocked_true_df" ]]; then
    echo "ERROR: Invalid -m value '$METHOD'. Use 'direct', 'cholesky', 'block_cholesky', 'true_df', or 'blocked_true_df'."
    exit 2
fi

ROSE_ORDER_LC="$(printf '%s' "$ROSE_ORDER" | tr '[:upper:]' '[:lower:]')"
if [[ "$ROSE_ORDER_LC" != "fragment" && "$ROSE_ORDER_LC" != "fragment_recanon" && "$ROSE_ORDER_LC" != "energy" ]]; then
    echo "ERROR: Invalid -o value '$ROSE_ORDER'. Use 'fragment', 'fragment_recanon', or 'energy'."
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

# Auto-cap default threads for large systems to avoid memory-pressure crashes
# (SIGBUS/Killed) during SCF/localization. Respect explicit -t from user.
if ! $USER_SET_OPENBLAS_THREADS || ! $USER_SET_OMP_THREADS; then
    if [[ -f "$DATA_DIR/mol_info.txt" ]]; then
        NAO_SUPER="$(sed -n '4p' "$DATA_DIR/mol_info.txt" | tr -d '[:space:]')"
        if [[ "$NAO_SUPER" =~ ^[0-9]+$ ]] && [[ "$NAO_SUPER" -ge 400 ]]; then
            echo "==> Large AO system detected (nao=$NAO_SUPER)."
            if ! $USER_SET_OPENBLAS_THREADS; then
                OPENBLAS_THREADS=1
            fi
            if ! $USER_SET_OMP_THREADS; then
                OMP_THREADS=1
            fi
            echo "    Auto-thread cap: OpenBLAS=$OPENBLAS_THREADS, OpenMP=$OMP_THREADS"
        fi
    fi
fi

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
echo "  LMO ordering : $ROSE_ORDER_LC"
echo "  OpenBLAS thd : $OPENBLAS_THREADS"
echo "  OpenMP thd   : $OMP_THREADS"
echo "  Data dir     : $DATA_DIR"
echo "------------------------------------------------------------"
EXAGRAD_FOCK_BUILDER="$METHOD_LC" \
EXAGRAD_ROSE_ORBITAL_ORDER="$ROSE_ORDER_LC" \
EXAGRAD_MOL_DIR="$DATA_DIR" \
OMP_NUM_THREADS="$OMP_THREADS" \
OPENBLAS_NUM_THREADS="$OPENBLAS_THREADS" \
./rhf_rose_main
echo "------------------------------------------------------------"

echo "==> Done."
