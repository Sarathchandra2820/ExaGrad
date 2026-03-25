import os
import numpy as np
from pyscf import gto, scf

def generate_integrals():
    # Define a simple molecule (Water in this case)
    # You can change the atoms and basis set as needed for your specific calculations.
    mol = gto.M(
        atom = '''
            O 0.0000000 0.0000000 0.0000000
            H 0.0000000 0.7570000 0.5870000
            H 0.0000000 -0.7570000 0.5870000
        ''',
        basis = 'sto-3g',
        verbose = 4,
        symmetry = False
    )

    # Run a self-consistent field (SCF) calculation
    # We use Restricted Hartree-Fock (RHF) as a starting point.
    mf = scf.RHF(mol)
    mf.kernel()

    # 1. 1-electron integrals (Core Hamiltonian = Kinetic + Nuclear Attraction)
    h_core = mol.intor('int1e_kin') + mol.intor('int1e_nuc')

    # 2. 2-electron integrals
    # Note: These are in the atomic orbital (AO) basis and are in the Chemists' notation (ij|kl).
    eri = mol.intor('int2e')

    # 3. Dipole integrals
    # Returns an array with shape (3, N_ao, N_ao) for x, y, z components.
    dipole_ints = mol.intor('int1e_r')

    # 4. MO coefficient matrix
    mo_coeff = mf.mo_coeff

    # --------------------------------------------------------------------------
    # Write the quantities to text files so Fortran can read them easily.
    # --------------------------------------------------------------------------
    
    os.makedirs('ints', exist_ok=True)
    
    # Save the Core Hamiltonian (2D array)
    np.savetxt('ints/hcore.txt', h_core, fmt='%24.16e')
    
    # Save the 2-electron integrals.
    # We flatten it to a 1D array to make it easier to read sequentially in Fortran.
    # The original shape is (N_ao, N_ao, N_ao, N_ao).
    np.savetxt('ints/eri.txt', eri.flatten(), fmt='%24.16e')
    
    # Save the dipole integrals (all components x, y, z)
    # We flatten it to a 1D array (shape 3 * N_ao * N_ao)
    np.savetxt('ints/dipole.txt', dipole_ints.flatten(), fmt='%24.16e')
    
    # Save the MO Coefficients (2D array)
    np.savetxt('ints/mo_coeff.txt', mo_coeff, fmt='%24.16e')

    print("\n--- Integral Generation Complete ---")
    print("Files written to the 'ints' directory:")
    print("- hcore.txt     (1-electron integrals)")
    print("- eri.txt       (2-electron integrals, flattened)")
    print("- dipole.txt    (dipole integrals, flattened)")
    print("- mo_coeff.txt  (MO coefficient matrix)")

if __name__ == '__main__':
    generate_integrals()
