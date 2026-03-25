#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

import numpy as np
from pyscf import gto, scf


def parse_fortran_mu(log_path: Path) -> np.ndarray:
    text = log_path.read_text()
    matches = re.findall(
        r"MU_ELEC_AU\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        text,
    )
    if not matches:
        raise ValueError("Could not find 'MU_ELEC_AU =' in Fortran output log")
    return np.array([float(matches[-1][0]), float(matches[-1][1]), float(matches[-1][2])])


def compute_pyscf_mu_elec_au(xyz: str, basis: str) -> np.ndarray:
    mol = gto.M(atom=xyz, basis=basis, cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    dm = mf.make_rdm1()
    dip_ao = mol.intor("int1e_r_sph", comp=3)
    mu = -np.einsum("xij,ji->x", dip_ao, dm)
    return mu


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare Fortran electronic dipole (a.u.) to PySCF")
    parser.add_argument("--xyz", required=True)
    parser.add_argument("--basis", required=True)
    parser.add_argument("--fortran-log", required=True)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()

    try:
        mu_f = parse_fortran_mu(Path(args.fortran_log))
    except Exception as exc:
        print(f"ERROR: {exc}")
        return 2

    try:
        mu_p = compute_pyscf_mu_elec_au(args.xyz, args.basis)
    except Exception as exc:
        print(f"ERROR: Failed to compute PySCF dipole: {exc}")
        return 2

    d = mu_f - mu_p
    ad = np.abs(d)

    print("==> Electronic dipole comparison (a.u.)")
    print(f"    Fortran : [{mu_f[0]: .10f}, {mu_f[1]: .10f}, {mu_f[2]: .10f}]")
    print(f"    PySCF   : [{mu_p[0]: .10f}, {mu_p[1]: .10f}, {mu_p[2]: .10f}]")
    print(f"    Delta   : [{d[0]: .3e}, {d[1]: .3e}, {d[2]: .3e}]")

    if np.max(ad) <= args.tolerance:
        print(f"    Status  : PASS (max |Delta| <= {args.tolerance:.1e})")
        return 0

    print(f"    Status  : WARN (max |Delta| > {args.tolerance:.1e})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
