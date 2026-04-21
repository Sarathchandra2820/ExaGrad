#!/usr/bin/env python3
"""
MolOrb — 3D Molecular Orbital Viewer for the Terminal
Real quantum chemistry (HF/3-21G) rendered as colored ASCII art.
"""

import curses
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto, scf

# ─── ASCII palette: low → high density ─────────────────────────────────────
CHARS = " ·:;+!?xX$#@"

# ─── Element → curses color-pair index ──────────────────────────────────────
ELEM_CP = {"C": 6, "N": 3, "O": 2, "S": 5, "F": 4, "P": 5, "Cl": 4, "Br": 2}

# ─── Common molecule name → SMILES ──────────────────────────────────────────
NAME_MAP = {
    "water": "O", "h2o": "O",
    "benzene": "c1ccccc1",
    "methane": "C", "ch4": "C",
    "ethanol": "CCO",
    "caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "co2": "O=C=O",
    "ammonia": "N", "nh3": "N",
    "formaldehyde": "C=O",
    "acetylene": "C#C",
    "ethylene": "C=C",
    "pyridine": "c1ccncc1",
    "phenol": "Oc1ccccc1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "naphthalene": "c1ccc2ccccc2c1",
}


def smiles_to_atoms(smiles: str) -> tuple:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result == -1:
        raise ValueError("Could not generate 3D coordinates.")
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
    return "\n".join(lines), mol


def run_hf(smiles: str):
    atom_str, rdmol = smiles_to_atoms(smiles)
    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = "3-21g"
    mol.verbose = 0
    mol.build()
    mf = scf.RHF(mol)
    mf.kernel()
    return mol, mf, rdmol


def eval_orbital(mol, mf, mo_idx: int, grid: int = 80, pad: float = 3.2):
    coords = mol.atom_coords()
    lo = coords.min(axis=0) - pad
    hi = coords.max(axis=0) + pad
    xs = np.linspace(lo[0], hi[0], grid)
    ys = np.linspace(lo[1], hi[1], grid)
    zs = np.linspace(lo[2], hi[2], grid)
    xx, yy, zz = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    ao = mol.eval_gto("GTOval", pts)
    mo_vals = ao @ mf.mo_coeff[:, mo_idx]
    return mo_vals.reshape(grid, grid, grid), lo, hi


def project_skeleton(mol, rdmol, R, W, H, grid, lo, hi):
    """Project heavy atoms and bonds into screen space using the same transform as render().
    Returns:
        atoms:     list of (sx, sy, label_char, color_pair_idx)
        bond_pts:  list of (sx, sy) for dots sampled along each bond
    """
    coords = mol.atom_coords()          # (natom, 3) Bohr, same order as rdmol atoms
    scale  = grid * 0.44
    aspect = 0.48
    center = grid / 2.0

    def _project(xyz_bohr):
        # Bohr → fractional grid index, centered, then apply rotation
        g  = (xyz_bohr - lo) / (hi - lo) * (grid - 1) - center
        p  = R @ (g / scale)
        return p[0] * scale + W / 2, p[1] * scale * aspect + H / 2

    atom_syms = [a.GetSymbol() for a in rdmol.GetAtoms()]
    heavy     = [i for i, s in enumerate(atom_syms) if s != "H"]
    heavy_set = set(heavy)

    atoms = []
    for i in heavy:
        sx, sy = _project(coords[i])
        sym = atom_syms[i]
        cp  = ELEM_CP.get(sym, 7)
        atoms.append((round(sx), round(sy), sym[0], cp))

    bond_pts = []
    for bond in rdmol.GetBonds():
        ai, aj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if ai not in heavy_set or aj not in heavy_set:
            continue
        ci, cj = coords[ai], coords[aj]
        # sample 12 intermediate points along the bond (endpoints drawn as atoms)
        for k in range(1, 12):
            t   = k / 12
            sx, sy = _project(ci * (1 - t) + cj * t)
            bond_pts.append((round(sx), round(sy)))

    return atoms, bond_pts


def rot_matrix(rx, ry, rz):
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def render(volume, R, W, H, iso=0.07):
    """Raycast orbital volume → 2D ASCII grid (vectorised over screen pixels)."""
    nx, ny, nz = volume.shape
    vmax = np.abs(volume).max()
    if vmax < 1e-10:
        return [[(0, " ")] * W for _ in range(H)]
    vol = volume / vmax
    scale = min(nx, ny, nz) * 0.44
    aspect = 0.48
    cx, cy, cz = nx / 2, ny / 2, nz / 2
    RT = R.T

    # Screen-space coordinates — shape (H, W)
    px = (np.arange(W) - W / 2) / scale
    py = (np.arange(H) - H / 2) / scale / aspect
    PX, PY = np.meshgrid(px, py)          # both (H, W)

    best      = np.zeros((H, W), dtype=np.float32)
    best_sign = np.zeros((H, W), dtype=np.int8)

    for ti in range(nz):
        t = (ti - nz / 2) / scale
        # Rotate each screen ray into volume space (no Python loop over pixels)
        ix = (RT[0, 0] * PX + RT[0, 1] * PY + RT[0, 2] * t) * scale + cx
        iy = (RT[1, 0] * PX + RT[1, 1] * PY + RT[1, 2] * t) * scale + cy
        iz = (RT[2, 0] * PX + RT[2, 1] * PY + RT[2, 2] * t) * scale + cz

        ixi = ix.astype(np.int32)
        iyi = iy.astype(np.int32)
        izi = iz.astype(np.int32)

        valid = (ixi >= 0) & (ixi < nx) & (iyi >= 0) & (iyi < ny) & (izi >= 0) & (izi < nz)
        ixi_c = np.clip(ixi, 0, nx - 1)
        iyi_c = np.clip(iyi, 0, ny - 1)
        izi_c = np.clip(izi, 0, nz - 1)

        v = np.where(valid, vol[ixi_c, iyi_c, izi_c], 0.0)
        improve = valid & (np.abs(v) > iso) & (np.abs(v) > np.abs(best))
        best      = np.where(improve, v, best)
        best_sign = np.where(improve, np.sign(v).astype(np.int8), best_sign)

    # Build result list
    screen = [[(0, " ")] * W for _ in range(H)]
    hit = np.abs(best) > iso
    if hit.any():
        char_idx = np.clip(
            ((np.abs(best) - iso) / (1.0 - iso) * (len(CHARS) - 1)).astype(int),
            0, len(CHARS) - 1,
        )
        ys, xs = np.where(hit)
        for sy, sx in zip(ys, xs):
            screen[sy][sx] = (int(best_sign[sy, sx]), CHARS[char_idx[sy, sx]])
    return screen


def draw_box(win, y, x, h, w, color):
    win.addstr(y,     x,     "╔" + "═" * (w - 2) + "╗", color)
    for row in range(1, h - 1):
        win.addstr(y + row, x, "║", color)
        win.addstr(y + row, x + w - 1, "║", color)
    win.addstr(y + h - 1, x, "╚" + "═" * (w - 2) + "╝", color)


def input_screen(stdscr):
    curses.curs_set(1)
    curses.echo()
    h, w = stdscr.getmaxyx()
    C = lambda n: curses.color_pair(n)

    stdscr.clear()
    title = "⚛   M O L O R B   ⚛"
    sub   = "3D Molecular Orbitals  ·  HF/3-21G  ·  ASCII Renderer"

    bx, by, bw, bh = w // 2 - 32, h // 2 - 7, 64, 14
    draw_box(stdscr, by, bx, bh, bw, C(1) | curses.A_BOLD)

    stdscr.addstr(by + 1, bx + (bw - len(title)) // 2, title,  C(5) | curses.A_BOLD)
    stdscr.addstr(by + 2, bx + (bw - len(sub)) // 2,   sub,    C(4))
    stdscr.addstr(by + 3, bx + 1, "─" * (bw - 2), C(1))

    examples = [
        ("water", "benzene", "caffeine", "aspirin"),
        ("ammonia", "ethylene", "pyridine", "naphthalene"),
    ]
    stdscr.addstr(by + 4, bx + 3, "Examples:", C(7) | curses.A_BOLD)
    for ri, row in enumerate(examples):
        line = "  ".join(f"[{m}]" for m in row)
        stdscr.addstr(by + 5 + ri, bx + 3, line, C(6))

    stdscr.addstr(by + 8,  bx + 3, "Enter molecule name or SMILES:", C(4) | curses.A_BOLD)
    stdscr.addstr(by + 9,  bx + 3, "❯ ", C(5) | curses.A_BOLD)
    stdscr.addstr(by + 11, bx + 3, "( ← / → ) rotate   ( ↑ / ↓ ) orbital   ( A ) spin", C(6))
    stdscr.addstr(by + 12, bx + 3, "( + / - ) isovalue  ( Q ) quit", C(6))
    stdscr.refresh()

    raw = stdscr.getstr(by + 9, bx + 5, 48).decode("utf-8").strip()
    curses.noecho()
    curses.curs_set(0)
    return raw


def loading_screen(stdscr, label):
    h, w = stdscr.getmaxyx()
    stdscr.clear()
    msg = f"  ⚛  {label}  "
    stdscr.addstr(h // 2 - 1, w // 2 - len(msg) // 2, msg,
                  curses.color_pair(5) | curses.A_BOLD)
    stdscr.addstr(h // 2 + 1, w // 2 - 18,
                  "  Please wait — running quantum chemistry...  ",
                  curses.color_pair(4))
    stdscr.refresh()


def orbital_label(idx, homo):
    d = idx - homo
    return {0: "HOMO", -1: "HOMO-1", -2: "HOMO-2",
            1: "LUMO", 2: "LUMO+1"}.get(d, f"MO {idx}")


def main(stdscr):
    curses.start_color()
    curses.use_default_colors()
    curses.init_pair(1, curses.COLOR_CYAN,    -1)  # borders
    curses.init_pair(2, curses.COLOR_RED,     -1)  # + lobe
    curses.init_pair(3, 33 if curses.COLORS >= 256 else curses.COLOR_BLUE, -1)  # − lobe
    curses.init_pair(4, curses.COLOR_GREEN,   -1)  # labels
    curses.init_pair(5, curses.COLOR_YELLOW,  -1)  # highlights
    curses.init_pair(6, curses.COLOR_WHITE,   -1)  # normal
    curses.init_pair(7, curses.COLOR_MAGENTA, -1)  # accent
    curses.curs_set(0)

    raw = input_screen(stdscr)
    smiles = NAME_MAP.get(raw.lower(), raw)

    loading_screen(stdscr, f"Computing {raw.upper()}")
    try:
        mol, mf, rdmol = run_hf(smiles)
    except Exception as e:
        h, w = stdscr.getmaxyx()
        stdscr.clear()
        stdscr.addstr(h // 2, w // 2 - 25,
                      f"  ✗ Error: {str(e)[:50]}  ",
                      curses.color_pair(2) | curses.A_BOLD)
        stdscr.addstr(h // 2 + 2, w // 2 - 12,
                      "Press any key to exit.", curses.color_pair(6))
        stdscr.nodelay(False)
        stdscr.getch()
        return

    homo = mol.nelectron // 2 - 1
    orb_idxs = list(range(max(0, homo - 2), min(len(mf.mo_energy), homo + 3)))

    loading_screen(stdscr, "Building orbital grids")
    volumes = {}
    grid_lo = grid_hi = grid_n = None
    for i in orb_idxs:
        vol, lo, hi = eval_orbital(mol, mf, i)
        volumes[i] = vol
        grid_lo, grid_hi, grid_n = lo, hi, vol.shape[0]

    # ── Main render loop ────────────────────────────────────────────────────
    rx, ry, rz = 0.3, 0.5, 0.0
    auto_spin = True
    iso = 0.08
    cur = orb_idxs.index(homo)
    stdscr.nodelay(True)

    while True:
        if auto_spin:
            ry += 0.03

        h, w = stdscr.getmaxyx()
        RH = h - 7          # render height
        RW = w - 2          # render width

        R   = rot_matrix(rx, ry, rz)
        scr = render(volumes[orb_idxs[cur]], R, RW, RH, iso)

        stdscr.erase()

        # ── skeleton (drawn first so orbital pixels overwrite where they overlap)
        sk_atoms, sk_bonds = project_skeleton(
            mol, rdmol, R, RW, RH, grid_n, grid_lo, grid_hi
        )
        for sx, sy, ch, cp in sk_atoms:
            if 0 <= sy < RH and 0 <= sx < RW:
                try:
                    stdscr.addstr(sy + 1, sx + 1, ch,
                                  curses.color_pair(cp) | curses.A_BOLD)
                except curses.error:
                    pass
        for sx, sy in sk_bonds:
            if 0 <= sy < RH and 0 <= sx < RW:
                try:
                    stdscr.addstr(sy + 1, sx + 1, "·", curses.color_pair(1))
                except curses.error:
                    pass

        # ── top bar
        lbl  = orbital_label(orb_idxs[cur], homo)
        ev   = mf.mo_energy[orb_idxs[cur]] * 27.211
        spin = "↻ spin" if auto_spin else "  paused"
        bar  = f" ⚛ MolOrb  │  {raw.upper()}  │  {lbl}  │  {ev:+.3f} eV  │  {spin} "
        pad  = "─" * w
        try:
            stdscr.addstr(0, 0, pad, curses.color_pair(1))
            stdscr.addstr(0, max(0, w // 2 - len(bar) // 2), bar,
                          curses.color_pair(1) | curses.A_BOLD)
        except curses.error:
            pass

        # ── orbital pixels
        for sy in range(min(RH, h - 7)):
            for sx in range(min(RW, w - 1)):
                sign, ch = scr[sy][sx]
                if ch != " ":
                    color = (curses.color_pair(2) if sign > 0
                             else curses.color_pair(3))
                    try:
                        stdscr.addstr(sy + 1, sx + 1, ch,
                                      color | curses.A_BOLD)
                    except curses.error:
                        pass

        # ── bottom panel
        bot = h - 6
        try:
            stdscr.addstr(bot, 0, "─" * w, curses.color_pair(1))

            stdscr.addstr(bot + 1, 2, "● +lobe", curses.color_pair(2) | curses.A_BOLD)
            stdscr.addstr(bot + 1, 12, "● −lobe", curses.color_pair(3) | curses.A_BOLD)
            stdscr.addstr(bot + 1, 22, f"iso={iso:.2f}", curses.color_pair(4))
            stdscr.addstr(bot + 1, 32,
                          f"atoms={rdmol.GetNumHeavyAtoms()}  "
                          f"e⁻={mol.nelectron}  basis=3-21G  method=RHF",
                          curses.color_pair(7))

            ctrl = "←/→ rotate   ↑/↓ orbital   A spin   +/- iso   Q quit"
            stdscr.addstr(bot + 2, max(0, w // 2 - len(ctrl) // 2),
                          ctrl, curses.color_pair(6))

            # orbital strip
            x = 2
            stdscr.addstr(bot + 3, x, "Orbitals:", curses.color_pair(4))
            x += 11
            for i, oi in enumerate(orb_idxs):
                tag = f" {orbital_label(oi, homo)} "
                if i == cur:
                    stdscr.addstr(bot + 3, x, tag,
                                  curses.color_pair(5) | curses.A_BOLD | curses.A_REVERSE)
                else:
                    stdscr.addstr(bot + 3, x, tag, curses.color_pair(6))
                x += len(tag) + 1

            stdscr.addstr(bot + 4, 0, "─" * w, curses.color_pair(1))
        except curses.error:
            pass

        stdscr.refresh()

        # ── input
        key = stdscr.getch()
        if key in (ord("q"), ord("Q")):
            break
        elif key == curses.KEY_LEFT:
            ry -= 0.15; auto_spin = False
        elif key == curses.KEY_RIGHT:
            ry += 0.15; auto_spin = False
        elif key == curses.KEY_UP:
            rx -= 0.15; auto_spin = False
        elif key == curses.KEY_DOWN:
            rx += 0.15; auto_spin = False
        elif key in (ord("a"), ord("A")):
            auto_spin = not auto_spin
        elif key in (ord("+"), ord("=")):
            iso = min(0.30, iso + 0.01)
        elif key == ord("-"):
            iso = max(0.01, iso - 0.01)
        elif key in (curses.KEY_PPAGE, ord("[")):
            cur = max(0, cur - 1)
        elif key in (curses.KEY_NPAGE, ord("]")):
            cur = min(len(orb_idxs) - 1, cur + 1)

        curses.napms(45)


if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        pass