"""从紧凑描述构造 molecule，供测试 fixture 与枚举器使用。"""

from __future__ import annotations

from .molecule import ATOM_SLOTS, Atom, Bond, Direction, EndKind, Layer, Molecule


def build(
    up: int,
    down: int,
    up_edges: list[tuple[int, int]],
    down_edges: list[tuple[int, int]],
    cross: list[tuple[int, int]],
) -> Molecule:
    """构造一个两层 molecule。

    Args:
        up: MU 的 atom 数，编号 0..up-1
        down: MD 的 atom 数，编号 up..up+down-1（外部用 0..down-1 表示，内部偏移）
        up_edges: MU 内部 bond，(parent, child) 用 MU 本地编号
        down_edges: MD 内部 bond，(parent, child) 用 MD 本地编号
        cross: 跨层 bond，(up_local, down_local)；方向恒为 UP 是 parent
               （Def 1.1(2)：MD 中无 atom 是 MU 中 atom 的 parent）

    余下的端点自动补为 free end，方向按该 atom 已有 bond 推断：
    朝上的空位记 TOP，朝下的记 BOTTOM。每个 atom 至多 2 parent / 2 child，
    因此 free end 的方向分配为：top 空位 = 2 - #parents，bottom 空位 = 2 - #children。
    """
    atoms: dict[int, Atom] = {}
    for i in range(up):
        atoms[i] = Atom(i, Layer.UP)
    for j in range(down):
        atoms[up + j] = Atom(up + j, Layer.DOWN)

    bonds: list[Bond] = []
    for p, c in up_edges:
        bonds.append(Bond(p, c))
    for p, c in down_edges:
        bonds.append(Bond(up + p, up + c))
    for u, d in cross:
        bonds.append(Bond(u, up + d))

    mol = Molecule(atoms, bonds)

    # 补 free end：每个 atom 有 2 个 top 槽、2 个 bottom 槽
    for aid, atom in atoms.items():
        n_par = len(mol.parents(aid))
        n_chi = len(mol.children(aid))
        if n_par > 2 or n_chi > 2:
            raise ValueError(f"atom {aid}: parents={n_par}, children={n_chi} (max 2)")
        for _ in range(2 - n_par):
            atom.ends.append((EndKind.FREE, Direction.TOP))
        for _ in range(2 - n_chi):
            atom.ends.append((EndKind.FREE, Direction.BOTTOM))
        total = len(mol.bonds_at(aid)) + len(atom.ends)
        if total != ATOM_SLOTS:
            raise ValueError(f"atom {aid}: slots={total} != {ATOM_SLOTS}")

    mol.check_invariants()
    return mol
