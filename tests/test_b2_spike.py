"""B2/X3C spike 的 corrected-W 回归。

这些有限实例不是 hardness 证明。它们锁定旧 gadget 在独立 `W` DP 下的真实
得分，并证明当前 YES/SCRAMBLE 对照没有 proposed exact-cover threshold gap。
"""

from __future__ import annotations

import pytest

from molecule_cut.builders import build
from molecule_cut.classify import W_of, is_toy1, n_2conn
from molecule_cut.fast_vd import fast_v_d
from molecule_cut.gadgets import (
    hubs_2conn,
    three_set_gadget,
    two_set_geometry,
    two_sets_half,
    two_sets_overlap_cross,
    two_sets_scramble,
    two_sets_yes,
)


@pytest.mark.parametrize(
    ("m", "expected_w", "expected_rho"),
    [(1, 1, 1), (2, 3, 3), (3, 4, 5), (4, 5, 7)],
)
def test_hubs_corrected_w(m, expected_w, expected_rho):
    mol = hubs_2conn(m)
    assert not is_toy1(mol)
    assert n_2conn(mol) == m
    mol.check_invariants()
    assert W_of(mol) == expected_w
    assert mol.circuit_rank() == expected_rho


def test_three_set_gadget_corrected_w():
    full = three_set_gadget()
    assert n_2conn(full) == 1
    full.check_invariants()
    assert W_of(full) == 3
    assert full.circuit_rank() == 2


def test_two_set_controls_corrected_w_have_no_yes_scramble_gap():
    """真实 W 下 YES 与 SCRAMBLE 相同，当前 gadget 没有 exact-cover gap。"""
    yes = two_sets_yes()
    half = two_sets_half()
    scr = two_sets_scramble()
    for mol in (yes, half, scr):
        mol.check_invariants()
    assert yes.circuit_rank() == scr.circuit_rank() == 5
    assert half.circuit_rank() == 2
    assert (W_of(yes), W_of(half), W_of(scr)) == (7, 6, 7)
    assert (fast_v_d(yes), fast_v_d(half), fast_v_d(scr)) == (-25, -28, -25)


def test_overlapping_sets_hit_slot_cap():
    """cap-1 元素上重叠 3-集无法同时嵌入（X3C 全实例归约障碍）。"""
    with pytest.raises(ValueError, match="parents"):
        two_set_geometry(two_sets_overlap_cross())


def test_path3_cross_pattern_changes_w():
    """固定两路径树，跨层模式可改变 W（H3：编码自由度存在）。"""
    ue = [(0, 1), (1, 2)]
    de = [(0, 1), (1, 2)]
    w_id = W_of(build(3, 3, ue, de, [(0, 0), (1, 1), (2, 2)]))
    w_sparse = W_of(build(3, 3, ue, de, [(2, 2)]))
    assert w_id == 2
    assert w_sparse == 1
