#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

from pyscf import gto, scf


def parse_fortran_energy(log_path: Path) -> float:
    text = log_path.read_text()
    matches = re.findall(r"E_total\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)", text)
    if not matches:
        raise ValueError("Could not find 'E_total =' in Fortran output log.")
    return float(matches[-1])


def compute_pyscf_energy(xyz: str, basis: str) -> float:
    mol = gto.M(atom=xyz, basis=basis, cart=False, verbose=0)
    mf = scf.RHF(mol)
    return float(mf.kernel())


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare ExaGrad SCF energy with PySCF RHF energy.")
    parser.add_argument("--xyz", required=True, help="Path to geometry .xyz file")
    parser.add_argument("--basis", required=True, help="Basis set name, e.g. sto-6g")
    parser.add_argument("--fortran-log", required=True, help="Path to captured Fortran output log")
    parser.add_argument("--tolerance", type=float, default=1.0e-6, help="Absolute tolerance for energy difference")
    args = parser.parse_args()

    log_path = Path(args.fortran_log)
    if not log_path.exists():
        print(f"ERROR: Fortran log not found: {log_path}")
        return 2

    try:
        e_fortran = parse_fortran_energy(log_path)
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 2

    try:
        e_pyscf = compute_pyscf_energy(args.xyz, args.basis)
    except Exception as exc:
        print(f"ERROR: Failed to compute PySCF energy: {exc}")
        return 2

    diff = e_fortran - e_pyscf
    adiff = abs(diff)

    print("==> Energy comparison (ExaGrad vs PySCF)")
    print(f"    Fortran RHF : {e_fortran: .10f} Ha")
    print(f"    PySCF RHF   : {e_pyscf: .10f} Ha")
    print(f"    Delta       : {diff: .10e} Ha")

    if adiff <= args.tolerance:
        print(f"    Status      : PASS (|Delta| <= {args.tolerance:.1e})")
        return 0

    print(f"    Status      : WARN (|Delta| > {args.tolerance:.1e})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
