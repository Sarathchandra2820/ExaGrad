import numpy as np
from pyscf import gto
import sys
import os

def export_cint_data(xyz_file, basis_name, output_dir):
    # Parse xyz to get atom coords and types
    mol = gto.M(
        atom=xyz_file,
        basis=basis_name,
        cart=False # We use spherical basis functions
    )
    
    # Extract libcint arrays
    atm = mol._atm
    bas = mol._bas
    env = mol._env

    # Save them as plain text so Fortran can read them sequentially
    np.savetxt(f"{output_dir}/atm.txt", atm.flatten(), fmt='%d')
    np.savetxt(f"{output_dir}/bas.txt", bas.flatten(), fmt='%d')
    np.savetxt(f"{output_dir}/env.txt", env, fmt='%24.16e')
    
    with open(f"{output_dir}/mol_info.txt", "w") as f:
        f.write(f"{atm.shape[0]}\n") # num atoms (natm)
        f.write(f"{bas.shape[0]}\n") # num shells (nbas)
        f.write(f"{env.shape[0]}\n") # env array size
        f.write(f"{mol.nao}\n") # num AO (spherical)

    print(f"Exported libcint env data to {output_dir}")
    print(f"Number of Atoms: {atm.shape[0]}")
    print(f"Number of Shells: {bas.shape[0]}")
    print(f"Total AO Basis Functions: {mol.nao}")

if __name__ == '__main__':
    xyz = sys.argv[1] if len(sys.argv) > 1 else 'geometry/H2O.xyz'
    basis = sys.argv[2] if len(sys.argv) > 2 else 'sto-3g'
    
    os.makedirs('ints', exist_ok=True)
    export_cint_data(xyz, basis, 'ints')
