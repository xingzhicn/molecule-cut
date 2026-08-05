"""Molecule 类标签与内在量便捷查询（Toy I / #2conn / W）。

供复杂度 spike 与测试使用；不改变 SPEC 语义。
"""

from __future__ import annotations

from .exhaustive import w
from .molecule import Layer, Molecule


def cross_degree_up(mol: Molecule) -> dict[int, int]:
    """每个 MU atom 的跨层 bond 数。"""
    counts: dict[int, int] = {u: 0 for u in mol.layer_atoms(Layer.UP)}
    for b in mol.cross_bonds():
        counts[b.parent] = counts.get(b.parent, 0) + 1
    return counts


def n_2conn(mol: Molecule) -> int:
    """#2conn = 恰有 2 条跨层 bond 的 MU atom 数（SPEC §1.5）。"""
    return sum(1 for c in cross_degree_up(mol).values() if c == 2)


def is_toy1(mol: Molecule) -> bool:
    """Toy I：每个 MU atom 恰 1 条跨层 bond（#2conn=0 且 X=MU）。"""
    deg = cross_degree_up(mol)
    return bool(deg) and all(c == 1 for c in deg.values())


def is_simplified_two_layer_shape(mol: Molecule) -> bool:
    """弱检查：两层非空、跨层方向 UP→DOWN、无初始 fixed end。

    不验证「每层是树」（调用方用 build 时已保证树形边集）；
    完整 Def 1.1 合法性仍以 build + check_invariants 为准。
    """
    up = mol.layer_atoms(Layer.UP)
    down = mol.layer_atoms(Layer.DOWN)
    if not up or not down:
        return False
    for atom in mol.atoms.values():
        if atom.n_fixed() != 0:
            return False
    for b in mol.cross_bonds():
        if mol.atoms[b.parent].layer is not Layer.UP:
            return False
        if mol.atoms[b.child].layer is not Layer.DOWN:
            return False
    return True


def W_of(mol: Molecule) -> int | None:
    """`W(M) = max_C #{33}(C)` 的精确计算。"""
    return w(mol)
