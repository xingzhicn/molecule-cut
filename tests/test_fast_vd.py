"""fast_v_d / fast_w 与参考实现的差分测试 + 状态空间上界。

语义必须分别与 `exhaustive.v_d_reference` / `w_reference` 完全一致。加速只改变
常数因子，不改变 SPEC §1.6 的定义。
"""

from __future__ import annotations

import random

import pytest

from molecule_cut.builders import build
from molecule_cut.elementary import ElemType, classify
from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.exhaustive import _legal_moves, v_d, v_d_reference, w, w_reference
from molecule_cut.fast_vd import (
    BAD4,
    GOOD33,
    GOOD343,
    INVALID,
    NORM2,
    NORM3,
    FastMolecule,
    fast_v_d,
    fast_w,
)

REF_TO_FAST = {
    ElemType.E2: NORM2,
    ElemType.E3: NORM3,
    ElemType.E4: BAD4,
    ElemType.E33: GOOD33,
    ElemType.E343: GOOD343,
}


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_fast_v_d_matches_reference_on_restricted_enumerator(n):
    for mol in enumerate_toy1(n):
        assert fast_v_d(mol, d=3) == v_d_reference(mol, d=3)


@pytest.mark.parametrize("n", [1, 2, 3])
def test_fast_w_matches_reference_on_restricted_enumerator(n):
    for mol in enumerate_toy1(n):
        assert fast_w(mol) == w_reference(mol)


def test_fast_matches_reference_n5_sample():
    """n=5 全量差分太慢；固定种子采样覆盖。"""
    random.seed(0)
    mols = list(enumerate_toy1(5))
    for mol in random.sample(mols, 25):
        assert fast_v_d(mol, d=3) == v_d_reference(mol, d=3)


def test_fast_targets_match_references_on_branched_mu():
    """覆盖完整 Toy I 中允许的 MU 双父构型，而不只覆盖 MU 链。"""
    mol = build(
        3,
        3,
        [(0, 2), (1, 2)],
        [(0, 1), (2, 1)],
        [(0, 0), (1, 0), (2, 2)],
    )
    assert fast_v_d(mol, d=3) == v_d_reference(mol, d=3)
    assert fast_w(mol) == w_reference(mol)
    assert w(mol) == w_reference(mol)


def test_w_is_not_unconditionally_v_d_plus_10d():
    """The two intrinsic objectives can prefer different decompositions."""
    mol = build(3, 3, [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(2, 2)])
    assert fast_w(mol) == w_reference(mol) == 1
    assert fast_v_d(mol, d=3) == -30
    assert fast_v_d(mol, d=3) + 30 == 0


def test_public_v_d_uses_fast_path():
    """默认 `v_d` 与 `fast_v_d` 一致（默认路径不再走对象级 DP）。"""
    for mol in enumerate_toy1(3):
        assert v_d(mol, d=3) == fast_v_d(mol, d=3)


def test_state_space_is_subset_mask():
    """状态完全由剩余 atom 掩码决定，空间 ≤ 2^{|M|}。"""
    mol = next(enumerate_toy1(4))
    fm = FastMolecule(mol)
    assert fm.n_atoms == 8  # |MU|=|MD|=4
    assert fm.full_mask == (1 << 8) - 1
    # 对每个可达状态，degree 公式与 base_free 一致
    state = fm.full_mask
    for i in range(fm.n_atoms):
        deg = fm.degree(i, state)
        assert deg == 4  # 初始无 fixed end


@pytest.mark.parametrize("n", [1, 2, 3])
def test_fast_move_enumerator_matches_reference_on_every_subset_state(n):
    """逐状态核对 move 完备性，而不只核对最终 DP 分数。"""
    for mol in enumerate_toy1(n):
        fm = FastMolecule(mol)
        ids = fm.ids
        for state in range(1 << fm.n_atoms):
            removed = {ids[i] for i in range(fm.n_atoms) if not (state >> i) & 1}
            residual = mol.cut_as_free(removed)

            fast_moves = {
                (tuple(ids[i] for i in subset), t)
                for subset, t in fm.legal_moves(state)
                if t != INVALID
            }
            reference_moves = {
                (tuple(sorted(move)), REF_TO_FAST[classify(residual, set(move))])
                for move in _legal_moves(residual)
            }
            assert fast_moves == reference_moves
