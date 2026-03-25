import pyscf
import h5py
import numpy as np 
import time

from pyscf import gto, scf, ao2mo

mol = gto.M(atom = 'geometry/ethane.xyz', basis='DZ')
#mf = scf.RHF(mol)
#mf.kernal()

eri_ao = mol.intor('int2e')
h1e_ao = mol.intor('int1e_ovlp')  # Overlap integrals
h1e_kin = mol.intor('int1e_kin')  # Kinetic energy integrals
h1e_nuc = mol.intor('int1e_nuc')  # Nuclear attraction integrals


with h5py.File('integrals.H5','w') as f:
    f.create_dataset('two_electron_integrals',data=eri_ao)
    f.create_dataset('overlap_integrals',data=h1e_ao)
    f.create_dataset('kinetic_energy_integrals',data=h1e_kin)
    f.create_dataset('nuclear_integrals',data=h1e_nuc)


def read_hdf5(filename, dataset_name):
    with h5py.File(filename, 'r') as file:
        # Access the dataset
        dataset = file[dataset_name]
        
        # Read the data into a numpy array
        data = np.array(dataset)
        
        # Print the shape and a preview of the data
        print(f"Shape of {dataset_name}: {data.shape}")
        print("Data preview:")
        print(data[:2, :2, :2, :2] if len(data.shape) == 4 else data[:5, :5])  # Adjust preview based on dimensions

    return data

# Example usage
filename = "integrals.H5"

# For a 2D dataset
data_2d = read_hdf5(filename, "overlap_integrals")




