"""枚举器同构去重的独立验证（PLAN V1.2）。

`enumerate.canonical_form` 用的是自写的 AHU 风格编码。此前只对 tie-breaking
枚举做过对拍，**分子层面的去重从未独立验证**。这里用一个「明显正确但低效」的
判据做对拍：直接枚举 MD 的全部标号置换，判断两个 molecule 是否同构。

同构的定义（与 canonical_form 的设计一致）：MU 是链，其编号由结构唯一确定
（顶端到最低端），故只允许置换 MD 的编号；bond 集合在该置换下相同即同构。
"""

from itertools import permutations

import pytest

from molecule_cut.enumerate import (
    canonical_form,
    enumerate_toy1,
    enumerate_toy1_general,
    rooted_binary_trees,
)
from molecule_cut.molecule import Layer

# OEIS A000111（Euler zigzag / 交错排列数）：每节点至多 2 child 的递增树计数
EULER_ZIGZAG = {1: 1, 2: 1, 3: 2, 4: 5, 5: 16, 6: 61, 7: 272}


@pytest.mark.parametrize("n", sorted(EULER_ZIGZAG))
def test_tree_generator_matches_oeis_a000111(n):
    """必要非充分检查：单树生成器的计数应等于 Euler zigzag 数。"""
    assert len(list(rooted_binary_trees(n))) == EULER_ZIGZAG[n]


def _bond_signature(mol, perm: dict[int, int]) -> tuple:
    """在 MD 编号置换 perm 下的 bond 集合签名。"""
    return tuple(
        sorted((perm.get(b.parent, b.parent), perm.get(b.child, b.child)) for b in mol.bonds)
    )


def _isomorphic_bruteforce(m1, m2) -> bool:
    """暴力判同构：枚举 MD 编号的全部置换。仅用于小规模对拍。"""
    up1, up2 = m1.layer_atoms(Layer.UP), m2.layer_atoms(Layer.UP)
    dn1, dn2 = m1.layer_atoms(Layer.DOWN), m2.layer_atoms(Layer.DOWN)
    if len(up1) != len(up2) or len(dn1) != len(dn2):
        return False
    if len(m1.bonds) != len(m2.bonds):
        return False

    target = _bond_signature(m2, {})
    for p in permutations(dn2):
        perm = dict(zip(dn1, p, strict=True))
        if _bond_signature(m1, perm) == target:
            return True
    return False


@pytest.mark.parametrize("n", [2, 3, 4])
def test_no_two_enumerated_molecules_are_isomorphic(n):
    """去重的正确性方向一：输出中不含互相同构的一对。"""
    mols = list(enumerate_toy1(n))
    for i in range(len(mols)):
        for j in range(i + 1, len(mols)):
            assert not _isomorphic_bruteforce(mols[i], mols[j]), (
                f"n={n}: molecules #{i} and #{j} are isomorphic but both were emitted"
            )


@pytest.mark.parametrize("n,nd", [(2, 2), (3, 2), (3, 3), (4, 3)])
def test_canonical_form_is_isomorphism_invariant(n, nd):
    """去重的正确性方向二：canonical_form 相同 ⟺ 暴力判据认为同构。

    这同时排除了「编码过粗（把不同构的合并）」与「编码过细（漏并同构的）」。
    """
    mols = list(enumerate_toy1_general(n, nd))
    for i in range(len(mols)):
        for j in range(i + 1, len(mols)):
            same_cf = canonical_form(mols[i], n, nd) == canonical_form(mols[j], n, nd)
            iso = _isomorphic_bruteforce(mols[i], mols[j])
            assert same_cf == iso, (
                f"(n={n},|MD|={nd}) #{i} vs #{j}: canonical_form equal={same_cf} "
                f"but bruteforce isomorphic={iso}"
            )


@pytest.mark.parametrize("n,nd", [(2, 2), (3, 3), (4, 3)])
def test_enumeration_is_complete_vs_raw_generation(n, nd):
    """完备性：原始生成（未去重）的每个 molecule 都同构于某个输出。"""
    from itertools import product

    from molecule_cut.builders import build

    emitted = list(enumerate_toy1_general(n, nd))
    up_edges = [(i, i + 1) for i in range(n - 1)]

    for par in rooted_binary_trees(nd):
        down_edges = [(par[i], i) for i in range(1, nd)]
        for assign in product(range(nd), repeat=n):
            try:
                raw = build(n, nd, up_edges, down_edges, [(u, assign[u]) for u in range(n)])
            except ValueError:
                continue
            assert any(_isomorphic_bruteforce(raw, e) for e in emitted), (
                f"(n={n},|MD|={nd}) raw molecule tree={par} assign={assign} "
                f"is isomorphic to no emitted molecule"
            )
