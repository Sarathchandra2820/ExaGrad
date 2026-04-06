#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

import numpy as np
from pyscf import gto, scf


def parse_fortran_alpha(log_path: Path) -> tuple[np.ndarray, float]:
    text = log_path.read_text()

    pat = re.compile(
        r"Static Dipole Polarizability Tensor \(a\.u\.\)\s*=*\n"
        r"\s*x\s*y\s*z\n"
        r"\s*x\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\n"
        r"\s*y\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\n"
        r"\s*z\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        re.MULTILINE,
    )

    m = pat.search(text)
    if not m:
        raise ValueError("Could not find static polarizability tensor in Fortran log.")

    alpha = np.array([float(x) for x in m.groups()], dtype=float).reshape(3, 3)

    iso_pat = re.compile(r"Isotropic alpha \(a\.u\.\)\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)")
    mi = iso_pat.search(text)
    iso = float(mi.group(1)) if mi else float(np.trace(alpha) / 3.0)

    return alpha, iso


def compute_pyscf_alpha_finite_field(
    xyz: str,
    basis: str,
    method: str,
    auxbasis: str,
    field_delta: float,
    conv_tol: float,
) -> np.ndarray:
    mol = gto.M(atom=xyz, basis=basis, cart=False, verbose=0)
    dip_ao = mol.intor_symmetric("int1e_r", comp=3)
    hcore0 = scf.hf.get_hcore(mol)

    def solve_mu(field_vec: np.ndarray) -> np.ndarray:
        mf = scf.RHF(mol)
        if method == "true_df":
            mf = mf.density_fit(auxbasis=auxbasis)
        mf.conv_tol = conv_tol
        mf.max_cycle = 100
        mf.get_hcore = lambda *args: hcore0 + np.einsum("x,xij->ij", field_vec, dip_ao)
        e = mf.kernel()
        if not mf.converged:
            raise RuntimeError(f"PySCF SCF did not converge under field, E={e}")
        dm = mf.make_rdm1()
        return -np.einsum("xij,ji->x", dip_ao, dm)

    alpha = np.zeros((3, 3), dtype=float)
    for j in range(3):
        f = np.zeros(3, dtype=float)
        f[j] = field_delta
        mu_p = solve_mu(f)
        mu_m = solve_mu(-f)
        alpha[:, j] = (mu_p - mu_m) / (2.0 * field_delta)

    return alpha


def fmt_vec(v: np.ndarray) -> str:
    return f"[{v[0]: .10f}, {v[1]: .10f}, {v[2]: .10f}]"


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare ExaGrad static polarizability (a.u.) with PySCF finite-field RHF.")
    parser.add_argument("--xyz", required=True, help="Path to geometry .xyz file")
    parser.add_argument("--basis", required=True, help="Basis set name")
    parser.add_argument("--fortran-log", required=True, help="Path to captured Fortran output log")
    parser.add_argument("--fortran-method", default="direct", help="Fortran mode: direct|cholesky|block_cholesky|true_df")
    parser.add_argument("--auxbasis", default="weigend", help="Auxiliary basis for true_df PySCF reference")
    parser.add_argument("--field-delta", type=float, default=1.0e-4, help="Finite electric-field step (a.u.)")
    parser.add_argument("--conv-tol", type=float, default=1.0e-11, help="PySCF RHF convergence tolerance")
    parser.add_argument("--tolerance", type=float, default=2.0e-3, help="Absolute tolerance for max tensor element difference")
    args = parser.parse_args()

    log_path = Path(args.fortran_log)
    if not log_path.exists():
        print(f"ERROR: Fortran log not found: {log_path}")
        return 2

    method = args.fortran_method.lower()

    try:
        alpha_f, iso_f = parse_fortran_alpha(log_path)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}")
        return 2

    try:
        alpha_p = compute_pyscf_alpha_finite_field(
            args.xyz,
            args.basis,
            method,
            args.auxbasis,
            args.field_delta,
            args.conv_tol,
        )
    except (RuntimeError, ValueError, OSError) as exc:
        print(f"ERROR: Failed to compute PySCF polarizability: {exc}")
        return 2

    iso_p = float(np.trace(alpha_p) / 3.0)
    delta = alpha_f - alpha_p
    max_abs = float(np.max(np.abs(delta)))

    print("==> Polarizability comparison (ExaGrad vs PySCF, a.u.)")
    print(f"    Fortran diag : {fmt_vec(np.diag(alpha_f))}")
    print(f"    PySCF diag   : {fmt_vec(np.diag(alpha_p))}")
    print(f"    Delta diag   : {fmt_vec(np.diag(delta))}")
    print(f"    Fortran iso  : {iso_f: .10f}")
    print(f"    PySCF iso    : {iso_p: .10f}")
    print(f"    Delta iso    : {iso_f - iso_p: .10e}")
    print(f"    Max |Delta|  : {max_abs: .10e}")

    if max_abs <= args.tolerance:
        print(f"    Status       : PASS (max |Delta| <= {args.tolerance:.1e})")
        return 0

    print(f"    Status       : WARN (max |Delta| > {args.tolerance:.1e})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
