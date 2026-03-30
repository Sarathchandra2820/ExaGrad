#!/usr/bin/env python3
import os
import sys

import numpy as np
from pyscf import gto, df
from pyscf.gto import mole


def export_df_metadata(xyz_file: str, basis_name: str, output_dir: str, auxbasis: str = "weigend") -> None:
    mol = gto.M(atom=xyz_file, basis=basis_name, cart=False, verbose=0)
    auxmol = df.addons.make_auxmol(mol, auxbasis=auxbasis)

    atm_all, bas_all, env_all = mole.conc_env(mol._atm, mol._bas, mol._env, auxmol._atm, auxmol._bas, auxmol._env)

    os.makedirs(output_dir, exist_ok=True)

    np.savetxt(os.path.join(output_dir, "df_meta_atm.txt"), atm_all.reshape(-1), fmt="%d")
    np.savetxt(os.path.join(output_dir, "df_meta_bas.txt"), bas_all.reshape(-1), fmt="%d")
    np.savetxt(os.path.join(output_dir, "df_meta_env.txt"), env_all.reshape(-1), fmt="%24.16e")

    with open(os.path.join(output_dir, "df_meta_info.txt"), "w", encoding="ascii") as f:
        f.write(f"{mol.natm}\n")
        f.write(f"{mol.nbas}\n")
        f.write(f"{auxmol.natm}\n")
        f.write(f"{auxmol.nbas}\n")
        f.write(f"{atm_all.shape[0]}\n")
        f.write(f"{bas_all.shape[0]}\n")
        f.write(f"{env_all.shape[0]}\n")
        f.write(f"{mol.nao}\n")
        f.write(f"{auxmol.nao}\n")

    print(f"Exported true DF metadata to {output_dir}")
    print(f"AO functions (mol): {mol.nao}")
    print(f"Aux functions      : {auxmol.nao}")
    print(f"Combined atoms     : {atm_all.shape[0]}")
    print(f"Combined shells    : {bas_all.shape[0]}")
    print(f"Aux basis          : {auxbasis}")


if __name__ == "__main__":
    xyz = sys.argv[1] if len(sys.argv) > 1 else "geometry/H2O.xyz"
    basis = sys.argv[2] if len(sys.argv) > 2 else "sto-3g"
    aux = sys.argv[3] if len(sys.argv) > 3 else "weigend"
    export_df_metadata(xyz, basis, "ints", aux)
