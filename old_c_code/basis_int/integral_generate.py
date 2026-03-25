import h5py # type: ignore
import numpy as np
import time
import sys

from pyscf import gto, scf, ao2mo

mol = gto.M(atom = '/Users/sarath/Documents/Research/other_projects/CCGF/copy_CCGF/geometry/H20.xyz', basis='sto-3g')

def libcint_input_gen(mol):
    bas = mol._bas
    env = mol._env
    coords = mol.atom_coords()
    return bas,env,coords

bas,env,coords = libcint_input_gen(mol)

with h5py.File('basis_info.H5','w') as f:
    f.create_dataset('bas',data=bas)
    f.create_dataset('env',data=env)
    f.create_dataset('coords',data=coords)