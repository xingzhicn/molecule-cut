"""B2/X3C spike 的构造守卫。

此前把 `v_d + 10d` 当成 W，故所有涉及 W 阈值的负结论均已冻结，等待以
真正的 W 重跑。这里只保留构造合法性与小实例的对象语义回归。
"""

from __future__ import annotations

import pytest

from molecule_cut.builders import build
from molecule_cut.classify import W_of, is_toy1, n_2conn
from molecule_cut.gadgets import (
    hubs_2conn,
    three_set_gadget,
    two_set_geometry,
    two_sets_half,
    two_sets_overlap_cross,
    two_sets_scramble,
    two_sets_yes,
)


def test_hubs_are_not_toy1_and_have_2conn():
    mol = hubs_2conn(3)
    assert not is_toy1(mol)
    assert n_2conn(mol) == 3
    mol.check_invariants()


def test_three_set_gadget_constructs():
    full = three_set_gadget()
    assert n_2conn(full) == 1
    full.check_invariants()


def test_two_set_controls_construct():
    """YES/half/scramble 是待重跑 W 的输入，不预先断言它们的阈值关系。"""
    yes = two_sets_yes()
    half = two_sets_half()
    scr = two_sets_scramble()
    for mol in (yes, half, scr):
        mol.check_invariants()
    assert yes.circuit_rank() == scr.circuit_rank() == 5
    assert half.circuit_rank() == 2


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
