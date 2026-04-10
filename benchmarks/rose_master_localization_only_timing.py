#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import statistics
import subprocess
import time
from pathlib import Path


def run_once(rose_exe: Path, workdir: Path, env: dict[str, str], log_path: Path) -> float:
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log_file:
        subprocess.run(
            [str(rose_exe)],
            cwd=str(workdir),
            env=env,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
        )
    end = time.perf_counter()
    return end - start


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Localization-only timing for rose-master using precomputed H5 inputs."
    )
    parser.add_argument(
        "--workdir",
        default="/Users/sarath/Documents/Research/ExaGrad/benchmarks/rose_master_pyridine_dimer",
        help="Directory containing INPUT_GENIBO, MOLECULE.h5, frag0.h5, frag1.h5",
    )
    parser.add_argument(
        "--rose-exe",
        default="/Users/sarath/Documents/Research/ExaGrad/rose-master/bin/rose.x",
        help="Path to rose-master executable",
    )
    parser.add_argument("--repeats", type=int, default=3, help="Number of timed repetitions")
    parser.add_argument("--threads", type=int, default=1, help="OMP/OpenBLAS thread count")
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    rose_exe = Path(args.rose_exe).resolve()

    required = [workdir / "INPUT_GENIBO", workdir / "MOLECULE.h5", workdir / "frag0.h5", workdir / "frag1.h5"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    if not rose_exe.exists():
        raise FileNotFoundError(f"ROSE executable not found: {rose_exe}")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["OPENBLAS_NUM_THREADS"] = str(args.threads)

    print("=== ROSE localization-only timing ===")
    print(f"workdir   : {workdir}")
    print(f"rose_exe  : {rose_exe}")
    print(f"repeats   : {args.repeats}")
    print(f"threads   : {args.threads}")
    print("includes  : rose.x localization workflow with precomputed H5 files")
    print("excludes  : PySCF SCF/H5 generation")
    print()

    timings: list[float] = []
    for idx in range(1, args.repeats + 1):
        log_path = workdir / f"localization_only_run_{idx}.log"
        dt = run_once(rose_exe=rose_exe, workdir=workdir, env=env, log_path=log_path)
        timings.append(dt)
        print(f"run {idx:02d}: {dt:8.3f} s   log={log_path.name}")

    mean = statistics.fmean(timings)
    std = statistics.pstdev(timings) if len(timings) > 1 else 0.0
    print()
    print(f"mean: {mean:.3f} s")
    print(f"std : {std:.3f} s")


if __name__ == "__main__":
    main()
