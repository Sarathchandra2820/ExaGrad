import os
import sys

import numpy as np
from pyscf import gto


def resolve_supersystem_xyz(root_dir, xyz_file):
    root_dir = os.path.abspath(root_dir)
    supra_xyz = os.path.join(root_dir, 'supra_mol.xyz')
    if os.path.isfile(supra_xyz):
        return supra_xyz
    return xyz_file


def find_fragment_xyz(fragment_dir):
    fragment_dir = os.path.abspath(fragment_dir)
    fragment_name = os.path.basename(fragment_dir)
    preferred_xyz = os.path.join(fragment_dir, f'{fragment_name}.xyz')
    if os.path.isfile(preferred_xyz):
        return preferred_xyz

    xyz_candidates = sorted(
        os.path.join(fragment_dir, entry)
        for entry in os.listdir(fragment_dir)
        if entry.lower().endswith('.xyz')
    )
    if len(xyz_candidates) == 1:
        return xyz_candidates[0]
    if len(xyz_candidates) > 1:
        raise ValueError(
            f"Multiple XYZ files found in {fragment_dir}. Expected one or {fragment_name}.xyz"
        )
    return None


def prepare_output_dirs(root_dir):
    root_dir = os.path.abspath(root_dir)
    ints_dir = os.path.join(root_dir, 'ints')
    os.makedirs(ints_dir, exist_ok=True)
    return root_dir, ints_dir


def find_fragment_xyz_files(root_dir):
    root_dir = os.path.abspath(root_dir)
    fragment_xyz_files = []

    for entry in sorted(os.listdir(root_dir)):
        entry_path = os.path.join(root_dir, entry)
        if entry.lower().startswith('frag_') and os.path.isdir(entry_path):
            fragment_xyz = find_fragment_xyz(entry_path)
            if fragment_xyz is not None:
                fragment_xyz_files.append(fragment_xyz)
            continue

        if os.path.isfile(entry_path) and entry.lower().endswith('.xyz') and entry.lower().startswith('frag_'):
            fragment_xyz_files.append(entry_path)

    return fragment_xyz_files


def export_fragment_int_dirs(root_dir, basis_name):
    root_dir = os.path.abspath(root_dir)
    fragment_xyz_files = find_fragment_xyz_files(root_dir)

    if not fragment_xyz_files:
        print(f"No fragment XYZ files found in {root_dir}")
        return []

    exported_fragment_dirs = []
    for fragment_xyz in fragment_xyz_files:
        fragment_name = os.path.splitext(os.path.basename(fragment_xyz))[0]
        fragment_root = os.path.join(root_dir, fragment_name)
        print(f"Exporting fragment integrals for {fragment_name} -> {fragment_root}")
        export_cint_data(fragment_xyz, basis_name, fragment_root)
        exported_fragment_dirs.append(fragment_root)

    return exported_fragment_dirs



def export_cint_data(xyz_file, basis_name, root_dir):
    root_dir, ints_dir = prepare_output_dirs(root_dir)

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
    np.savetxt(os.path.join(ints_dir, 'atm.txt'), atm.flatten(), fmt='%d')
    np.savetxt(os.path.join(ints_dir, 'bas.txt'), bas.flatten(), fmt='%d')
    np.savetxt(os.path.join(ints_dir, 'env.txt'), env, fmt='%24.16e')

    with open(os.path.join(root_dir, 'mol_info.txt'), 'w', encoding='ascii') as f:
        f.write(f"{atm.shape[0]}\n") # num atoms (natm)
        f.write(f"{bas.shape[0]}\n") # num shells (nbas)
        f.write(f"{env.shape[0]}\n") # env array size
        f.write(f"{mol.nao}\n") # num AO (spherical)
        f.write(f"{mol.nelectron}\n") # num electrons

    print(f"Exported libcint env data to {root_dir}")
    print(f"Number of Atoms: {atm.shape[0]}")
    print(f"Number of Shells: {bas.shape[0]}")
    print(f"Total AO Basis Functions: {mol.nao}")

    # -----------------------------------------------------------------------
    # Export minimal basis (minao) overlap data for IAO/IBO localization.
    # S21 = <minao | full_basis> cross-overlap (nmin x nfull)
    # S22 = <minao | minao> self-overlap (nmin x nmin)
    # nmo_frag = number of minao functions per atom (for IBO fragment assignment)
    # -----------------------------------------------------------------------
    # mol_min = gto.M(
    #     atom=xyz_file,
    #     basis='minao',
    #     cart=False
    # )

    # S21 = gto.intor_cross('int1e_ovlp', mol_min, mol)  # (nmin, nfull)
    # S22 = mol_min.intor('int1e_ovlp')                  # (nmin, nmin)

    # # Count minao spherical AOs per atom
    # nmo_frag = []
    # for iatm in range(mol.natm):
    #     count = 0
    #     for ish in range(mol_min.nbas):
    #         if mol_min.bas_atom(ish) == iatm:
    #             l = mol_min.bas_angular(ish)
    #             count += 2 * l + 1  # spherical: 2l+1 per shell
    #     nmo_frag.append(count)

    # np.savetxt(os.path.join(ints_dir, 'S21.txt'), S21.flatten(), fmt='%24.16e')
    # np.savetxt(os.path.join(ints_dir, 'S22.txt'), S22.flatten(), fmt='%24.16e')

    # with open(os.path.join(ints_dir, 'minbas_info.txt'), 'w', encoding='ascii') as f:
    #     f.write(f"{mol_min.nao}\n")       # nmin: total minao functions
    #     f.write(f"{mol.natm}\n")           # natm: number of atoms
    #     for n in nmo_frag:
    #         f.write(f"{n}\n")             # nmo_frag(i) for each atom

    # print(f"Exported minao basis overlap data (nmin={mol_min.nao}, nmo_frag={nmo_frag})")

    

def write_rose_info(root_dir, fragment_dirs):
    """Write rose_info.txt manifest listing fragment count and directory names."""
    root_dir = os.path.abspath(root_dir)
    info_path = os.path.join(root_dir, 'rose_info.txt')
    with open(info_path, 'w', encoding='ascii') as f:
        f.write(f"{len(fragment_dirs)}\n")
        for frag_dir in fragment_dirs:
            f.write(f"{os.path.basename(frag_dir)}\n")
    print(f"Wrote rose manifest: {info_path} ({len(fragment_dirs)} fragments)")


if __name__ == '__main__':
    xyz = sys.argv[1] if len(sys.argv) > 1 else 'geometry/H2O.xyz'
    basis = sys.argv[2] if len(sys.argv) > 2 else 'sto-3g'
    output_root = sys.argv[3] if len(sys.argv) > 3 else '.'
    xyz = resolve_supersystem_xyz(output_root, xyz)

    export_cint_data(xyz, basis, output_root)
    frag_dirs = export_fragment_int_dirs(output_root, basis)
    if frag_dirs:
        write_rose_info(output_root, frag_dirs)
