#!/usr/bin/env python3
import os
import sys
import numpy as np
from pyscf import gto, df

from export_cint_env import find_fragment_xyz_files


def export_true_df_factors(xyz_file: str, basis_name: str, root_dir: str, auxbasis: str = "weigend") -> None:
    mol = gto.M(atom=xyz_file, basis=basis_name, cart=False, verbose=0)
    auxmol = df.addons.make_auxmol(mol, auxbasis=auxbasis)

    # cderi shape: (naux, nao*(nao+1)/2) where each row is packed lower-triangular AO pair.
    cderi = np.asarray(df.incore.cholesky_eri(mol, auxbasis=auxbasis), dtype=np.float64)
    nao = mol.nao
    naux = cderi.shape[0]
    npair = cderi.shape[1]

    root_dir = os.path.abspath(root_dir)
    output_dir = os.path.join(root_dir, "ints")
    os.makedirs(output_dir, exist_ok=True)

    # Lightweight hybrid export: aux environment + packed cderi.
    np.savetxt(os.path.join(output_dir, "aux_atm.txt"), auxmol._atm.flatten(), fmt="%d")
    np.savetxt(os.path.join(output_dir, "aux_bas.txt"), auxmol._bas.flatten(), fmt="%d")
    np.savetxt(os.path.join(output_dir, "aux_env.txt"), auxmol._env, fmt="%24.16e")

    with open(os.path.join(output_dir, "aux_info.txt"), "w", encoding="ascii") as f:
        f.write(f"{auxmol._atm.shape[0]}\n")
        f.write(f"{auxmol._bas.shape[0]}\n")
        f.write(f"{auxmol._env.shape[0]}\n")
        f.write(f"{auxmol.nao}\n")

    with open(os.path.join(output_dir, "df_info.txt"), "w", encoding="ascii") as f:
        f.write(f"{nao}\n")
        f.write(f"{naux}\n")
        f.write(f"{npair}\n")

    # Packed rows written in Q-major order, then pair index within each row.
    np.savetxt(os.path.join(output_dir, "df_cderi.txt"), cderi.reshape(-1), fmt="%24.16e")

    print(f"Exported true DF auxiliary data to {output_dir}")
    print(f"AO functions: {nao}")
    print(f"Aux functions: {naux}")
    print(f"AO-pairs (packed): {npair}")
    print(f"Aux basis: {auxbasis}")


def export_fragment_true_df_factors(root_dir: str, basis_name: str, auxbasis: str = "weigend") -> list[str]:
    root_dir = os.path.abspath(root_dir)
    fragment_xyz_files = find_fragment_xyz_files(root_dir)

    if not fragment_xyz_files:
        print(f"No fragment XYZ files found in {root_dir}")
        return []

    exported_fragment_dirs = []
    for fragment_xyz in fragment_xyz_files:
        fragment_name = os.path.splitext(os.path.basename(fragment_xyz))[0]
        fragment_root = os.path.join(root_dir, fragment_name)
        print(f"Exporting DF factors for {fragment_name} -> {fragment_root}")
        export_true_df_factors(fragment_xyz, basis_name, fragment_root, auxbasis)
        exported_fragment_dirs.append(fragment_root)

    return exported_fragment_dirs


if __name__ == "__main__":
    xyz = sys.argv[1] if len(sys.argv) > 1 else "geometry/H2O.xyz"
    basis = sys.argv[2] if len(sys.argv) > 2 else "sto-3g"
    aux = sys.argv[3] if len(sys.argv) > 3 else "weigend"
    output_root = sys.argv[4] if len(sys.argv) > 4 else "."
    export_true_df_factors(xyz, basis, output_root, aux)
    export_fragment_true_df_factors(output_root, basis, aux)
