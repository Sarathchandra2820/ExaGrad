#!/bin/bash
# ============================================================
#  ExaGrad — Compile & Run ROSE-RHF (supersystem + fragments)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON_CMD=("${PYTHON:-python3}")
CONDA_ENV="${CONDA_ENV:-exagrad}"
if ! "${PYTHON_CMD[@]}" -c "import numpy, pyscf" >/dev/null 2>&1; then
    if command -v conda >/dev/null 2>&1; then
        PYTHON_CMD=(conda run -n "$CONDA_ENV" python)
    fi
fi

if ! "${PYTHON_CMD[@]}" -c "import numpy, pyscf" >/dev/null 2>&1; then
    echo "ERROR: Could not find a Python interpreter with numpy and pyscf."
    echo "Set PYTHON to a working interpreter or set CONDA_ENV to the correct conda environment."
    exit 2
fi

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
USER_SET_METHOD=false
NVAL_VIR_LIST=""
USER_SET_NVAL_VIR=false

# --- parse optional arguments --------------------------------
usage() {
    echo "Usage: $0 [-g geometry.xyz] [-b basis_set] [-m direct|cholesky|block_cholesky|true_df|blocked_true_df] [-a auxbasis] [-d data_dir] [-t nthreads] [-t1 openblas_threads] [-t2 openmp_threads] [-v nval_list] [-o fragment|fragment_recanon|energy] [-c] [-h]"
    echo ""
    echo "  -g FILE    Path to supersystem .xyz file     (default: DATA_DIR/supra_mol.xyz if present)"
    echo "  -b BASIS   Basis set name                    (default: $BASIS)"
    echo "  -m METHOD  Fock builder: direct|cholesky|block_cholesky|true_df|blocked_true_df  (default: $METHOD)"
    echo "  -a AUX     Auxiliary basis for true_df       (default: $AUXBASIS)"
    echo "  -d DIR     Molecule data root directory      (default: $DATA_DIR)"
    echo "  -t N       Set both OpenBLAS and OpenMP threads to N (default: $DEFAULT_THREADS)"
    echo "  -t1 N      Set OpenBLAS threads only         (default: $DEFAULT_THREADS)"
    echo "  -t2 N      Set OpenMP threads only           (default: $DEFAULT_THREADS)"
    echo "  -v LIST    Comma-separated n_val_vir per fragment in rose_info order"
    echo "             Example: -v 3,4,2 (use -1 for 'all virtuals')"
    echo "  -o MODE    LMO ordering: fragment|fragment_recanon|energy (default: fragment)"
    echo "  -c         Clean build (make clean first)"
    echo "  -h         Show this help message"
    exit 0
}

apply_nval_vir_overrides() {
    local data_dir="$1"
    local rose_info="$data_dir/rose_info.txt"
    local nfrags frag_line frag_name val idx
    local -a rose_lines
    local -a nval_values

    if ! $USER_SET_NVAL_VIR; then
        return 0
    fi

    if [[ ! -f "$rose_info" ]]; then
        echo "ERROR: rose_info.txt not found at '$rose_info'."
        exit 2
    fi

    rose_lines=()
    while IFS= read -r frag_line; do
        rose_lines+=("$frag_line")
    done < "$rose_info"
    if [[ ${#rose_lines[@]} -lt 1 ]]; then
        echo "ERROR: rose_info.txt is empty."
        exit 2
    fi

    nfrags="$(printf '%s' "${rose_lines[0]}" | tr -d '[:space:]')"
    if ! [[ "$nfrags" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Invalid fragment count in rose_info.txt: '${rose_lines[0]}'."
        exit 2
    fi

    IFS=',' read -r -a nval_values <<< "$NVAL_VIR_LIST"
    if [[ ${#nval_values[@]} -ne "$nfrags" ]]; then
        echo "ERROR: -v expects exactly $nfrags values, got ${#nval_values[@]}."
        exit 2
    fi

    for ((idx=0; idx<nfrags; idx++)); do
        val="$(printf '%s' "${nval_values[$idx]}" | tr -d '[:space:]')"
        if ! [[ "$val" =~ ^-?[0-9]+$ ]]; then
            echo "ERROR: Invalid -v entry '${nval_values[$idx]}' (must be integer >= -1)."
            exit 2
        fi
        if (( val < -1 )); then
            echo "ERROR: Invalid -v entry '$val' (minimum allowed is -1)."
            exit 2
        fi
        nval_values[$idx]="$val"
    done

    {
        echo "$nfrags"
        for ((idx=1; idx<=nfrags; idx++)); do
            if [[ $idx -ge ${#rose_lines[@]} ]]; then
                echo "ERROR: rose_info.txt has fewer fragment lines than expected ($nfrags)." >&2
                exit 2
            fi
            frag_line="${rose_lines[$idx]}"
            frag_name="$(printf '%s\n' "$frag_line" | awk '{print $1}')"
            if [[ -z "$frag_name" ]]; then
                echo "ERROR: Could not parse fragment name from rose_info line $((idx+1))." >&2
                exit 2
            fi
            echo "$frag_name ${nval_values[$((idx-1))]}"
        done
    } > "$rose_info"

    echo "==> Applied n_val_vir overrides from -v to $rose_info"
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
            METHOD="$2"; USER_SET_METHOD=true; shift 2 ;;
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
        -v)
            NVAL_VIR_LIST="$2"
            USER_SET_NVAL_VIR=true
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
"${PYTHON_CMD[@]}" python/export_cint_env.py "$XYZ" "$BASIS" "$DATA_DIR"
apply_nval_vir_overrides "$DATA_DIR"

# Auto-cap resources for large systems to avoid memory-pressure crashes
# (SIGBUS/Killed) during SCF/localization. Thread caps respect explicit -t flags;
# method auto-switch respects explicit -m flag.
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
        # true_df stores the full dense B(nao,nao,naux) tensor; for nao>=400 this
        # can require >7 GB peak (rhs + B_dense simultaneously) and segfaults.
        # blocked_true_df stores only the triangular-packed Bp(naux,npair) tensor
        # (~half the peak, ~half the live footprint) and contracts J/K in blocks.
        if [[ "$METHOD_LC" == "true_df" ]] && ! $USER_SET_METHOD; then
            METHOD_LC="blocked_true_df"
            echo "    Auto-switched method: true_df -> blocked_true_df"
            echo "      (dense CDERI tensor ~${NAO_SUPER}^2 * naux * 8 bytes would exhaust memory)"
        fi
    fi
fi

if [[ "$METHOD_LC" == "true_df" ]]; then
    echo "==> Exporting true DF metadata (supersystem + fragments) …"
    "${PYTHON_CMD[@]}" python/export_df_metadata.py "$XYZ" "$BASIS" "$AUXBASIS" "$DATA_DIR"
elif [[ "$METHOD_LC" == "blocked_true_df" ]]; then
    echo "==> Exporting blocked true-DF packed factors (supersystem + fragments) …"
    "${PYTHON_CMD[@]}" python/export_df_metadata.py "$XYZ" "$BASIS" "$AUXBASIS" "$DATA_DIR"
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
