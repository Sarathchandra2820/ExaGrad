#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any


def parse_xyz(xyz_path: Path) -> list[tuple[str, tuple[float, float, float]]]:
    lines = xyz_path.read_text(encoding="utf-8").strip().splitlines()
    natoms = int(lines[0].strip())
    atom_lines = lines[2:2 + natoms]
    atoms: list[tuple[str, tuple[float, float, float]]] = []
    for line in atom_lines:
        parts = line.split()
        atoms.append((parts[0], (float(parts[1]), float(parts[2]), float(parts[3]))))
    return atoms


def split_fragments(atoms: list[tuple[str, tuple[float, float, float]]], nfrags: int, natoms_per_frag: int | None) -> list[list[tuple[str, tuple[float, float, float]]]]:
    if natoms_per_frag is None:
        if len(atoms) % nfrags != 0:
            raise ValueError("Cannot split atoms evenly; provide --natoms-per-frag")
        natoms_per_frag = len(atoms) // nfrags
    if natoms_per_frag * nfrags > len(atoms):
        raise ValueError("Requested fragments exceed total atoms")

    fragments: list[list[tuple[str, tuple[float, float, float]]]] = []
    for idx in range(nfrags):
        start = idx * natoms_per_frag
        end = start + natoms_per_frag
        fragments.append(atoms[start:end])
    return fragments


def write_input_genibo(path: Path, nfrags: int, version: str, charge: int, exponent: int) -> None:
    content = (
        "**ROSE\n"
        ".VERSION\n"
        f"{version}\n"
        ".CHARGE\n"
        f"{charge}\n"
        ".EXPONENT\n"
        f"{exponent}\n"
        ".FILE_FORMAT\n"
        "h5\n"
        ".NFRAGMENTS\n"
        f"{nfrags}\n"
        "*END OF INPUT\n"
    )
    path.write_text(content, encoding="utf-8")


def format_seconds(seconds: float) -> str:
    return f"{seconds:9.3f} s"


def move_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        if dst.exists():
            dst.unlink()
        shutil.move(str(src), str(dst))


def main() -> None:
    parser = argparse.ArgumentParser(description="Full PySCF+ROSE pipeline timing with phase-by-phase wall-clock breakdown.")
    parser.add_argument("--xyz", default="/Users/sarath/Documents/Research/ExaGrad/geometry/pyridine_dimer.xyz", help="Supersystem XYZ path")
    parser.add_argument("--basis", default="cc-pVDZ", help="Basis set for supermolecule and fragments")
    parser.add_argument("--workdir", default="/Users/sarath/Documents/Research/ExaGrad/benchmarks/rose_master_pyridine_dimer_recalc", help="Benchmark work directory")
    parser.add_argument("--rose-root", default="/Users/sarath/Documents/Research/ExaGrad/rose-master", help="rose-master root path")
    parser.add_argument("--version", default="Stndrd_2013", help="ROSE .VERSION for INPUT_GENIBO")
    parser.add_argument("--exponent", type=int, default=2, help="ROSE localization exponent")
    parser.add_argument("--charge", type=int, default=0, help="Molecular charge")
    parser.add_argument("--multiplicity", type=int, default=1, help="Molecular multiplicity")
    parser.add_argument("--nfrags", type=int, default=2, help="Number of fragments")
    parser.add_argument("--natoms-per-frag", type=int, default=11, help="Atoms per fragment; set 0 to auto")
    parser.add_argument("--threads", type=int, default=1, help="OMP/OpenBLAS thread count")
    parser.add_argument("--report", default="timing_breakdown.json", help="Output JSON filename inside workdir")
    parser.add_argument("--keep-existing", action="store_true", help="Do not delete existing workdir files before run")
    args = parser.parse_args()

    xyz_path = Path(args.xyz).resolve()
    rose_root = Path(args.rose_root).resolve()
    rose_exe = rose_root / "bin" / "rose.x"
    workdir = Path(args.workdir).resolve()

    if not xyz_path.exists():
        raise FileNotFoundError(f"XYZ not found: {xyz_path}")
    if not rose_exe.exists():
        raise FileNotFoundError(f"ROSE executable not found: {rose_exe}")

    if workdir.exists() and not args.keep_existing:
        for item in workdir.iterdir():
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
    workdir.mkdir(parents=True, exist_ok=True)

    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)

    import sys
    sys.path.append(str(rose_root / "python_scripts" / "genibo"))
    from genibo_MOcoeff import moldata, run_pyscf_mocoeff

    atoms = parse_xyz(xyz_path)
    natoms_per_frag = None if args.natoms_per_frag == 0 else args.natoms_per_frag
    fragments = split_fragments(atoms, args.nfrags, natoms_per_frag)

    print("=== Full PySCF + ROSE pipeline timing ===")
    print(f"xyz        : {xyz_path}")
    print(f"basis      : {args.basis}")
    print(f"workdir    : {workdir}")
    print(f"threads    : {args.threads}")
    print(f"fragments  : {args.nfrags} (atoms/frag={len(fragments[0])})")
    print(f"rose exe   : {rose_exe}")
    print()

    timings: dict[str, float] = {}
    phase_meta: dict[str, Any] = {}

    t_total0 = time.perf_counter()

    t0 = time.perf_counter()
    write_input_genibo(workdir / "INPUT_GENIBO", args.nfrags, args.version, args.charge, args.exponent)
    timings["input_preparation"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    molecule = moldata(geometry=atoms, multiplicity=args.multiplicity, charge=args.charge, data_directory=str(workdir))
    molecule, pyscf_mol, mf = run_pyscf_mocoeff(
        molecule,
        basis=args.basis,
        relativistic=False,
        spherical=False,
        restricted=(args.multiplicity == 1),
        openshell=(args.multiplicity > 1),
        save=True,
        uncontract=False,
    )
    move_if_exists(workdir / f"{molecule.name}.h5", workdir / "MOLECULE.h5")
    move_if_exists(workdir / f"{molecule.name}.chk", workdir / "MOLECULE.chk")
    timings["pyscf_super_scf"] = time.perf_counter() - t0
    phase_meta["super_nao"] = int(pyscf_mol.nao_nr())
    phase_meta["super_nelec"] = int(pyscf_mol.nelectron)

    fragment_times: list[float] = []
    for idx, frag_geom in enumerate(fragments):
        t0 = time.perf_counter()
        fragment = moldata(geometry=frag_geom, multiplicity=1, charge=0, data_directory=str(workdir))
        fragment, frag_mol, _ = run_pyscf_mocoeff(
            fragment,
            basis=args.basis,
            relativistic=False,
            spherical=False,
            restricted=True,
            openshell=False,
            save=True,
            uncontract=False,
        )
        move_if_exists(workdir / f"{fragment.name}.h5", workdir / f"frag{idx}.h5")
        move_if_exists(workdir / f"{fragment.name}.chk", workdir / f"frag{idx}.chk")
        dt = time.perf_counter() - t0
        fragment_times.append(dt)
        phase_meta[f"frag{idx}_nao"] = int(frag_mol.nao_nr())
        phase_meta[f"frag{idx}_nelec"] = int(frag_mol.nelectron)

    timings["pyscf_fragments_total"] = sum(fragment_times)
    for idx, dt in enumerate(fragment_times):
        timings[f"pyscf_fragment_{idx}"] = dt

    t0 = time.perf_counter()
    rose_log = workdir / "rose_localization.log"
    with rose_log.open("w", encoding="utf-8") as out:
        subprocess.run(
            [str(rose_exe)],
            cwd=str(workdir),
            env=os.environ.copy(),
            stdout=out,
            stderr=subprocess.STDOUT,
            check=True,
        )
    timings["rose_localization"] = time.perf_counter() - t0

    timings["total_pipeline"] = time.perf_counter() - t_total0

    known_parts = (
        timings["input_preparation"]
        + timings["pyscf_super_scf"]
        + timings["pyscf_fragments_total"]
        + timings["rose_localization"]
    )
    timings["other_overhead"] = max(timings["total_pipeline"] - known_parts, 0.0)

    report = {
        "settings": {
            "xyz": str(xyz_path),
            "basis": args.basis,
            "workdir": str(workdir),
            "rose_executable": str(rose_exe),
            "version": args.version,
            "exponent": args.exponent,
            "charge": args.charge,
            "multiplicity": args.multiplicity,
            "nfrags": args.nfrags,
            "atoms_per_fragment": len(fragments[0]),
            "threads": args.threads,
        },
        "timings_seconds": timings,
        "metadata": phase_meta,
    }

    report_path = workdir / args.report
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("--- Timing breakdown (wall time) ---")
    print(f"input_preparation       : {format_seconds(timings['input_preparation'])}")
    print(f"pyscf_super_scf         : {format_seconds(timings['pyscf_super_scf'])}")
    for idx in range(args.nfrags):
        print(f"pyscf_fragment_{idx}        : {format_seconds(timings[f'pyscf_fragment_{idx}'])}")
    print(f"pyscf_fragments_total   : {format_seconds(timings['pyscf_fragments_total'])}")
    print(f"rose_localization       : {format_seconds(timings['rose_localization'])}")
    print(f"other_overhead          : {format_seconds(timings['other_overhead'])}")
    print(f"total_pipeline          : {format_seconds(timings['total_pipeline'])}")
    print()
    print(f"Report written to: {report_path}")
    print(f"ROSE log written to: {workdir / 'rose_localization.log'}")


if __name__ == "__main__":
    main()
