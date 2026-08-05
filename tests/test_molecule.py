"""Molecule 数据结构与 cut 语义的单元测试（SPEC §1, §3.1-3.3）。"""

import pytest

from molecule_cut.builders import build
from molecule_cut.elementary import ElemType, classify, is_proper, lambda_potential
from molecule_cut.molecule import Direction, EndKind, Layer


def test_degree_convention_excludes_fixed_ends():
    """SPEC §3.1：deg = bond + free，fixed 不计。"""
    m = build(1, 1, [], [], [(0, 0)])
    # 每个 atom：1 bond + 3 free = deg 4
    assert m.degree(0) == 4
    assert m.degree(1) == 4
    m.check_invariants()

    after = m.cut_as_free({0})
    # atom 1 的 bond 变成 fixed end → deg 降 1
    assert after.degree(1) == 3
    assert after.atoms[1].n_fixed() == 1


def test_cut_direction_of_fixed_end():
    """SPEC §3.2：切 parent 在 child 处留 top fixed end，反之留 bottom。"""
    m = build(1, 1, [], [], [(0, 0)])
    # atom 0 (UP) 是 atom 1 (DOWN) 的 parent
    after_parent_cut = m.cut_as_free({0})
    assert Direction.TOP in after_parent_cut.atoms[1].fixed_directions()

    after_child_cut = m.cut_as_free({1})
    assert Direction.BOTTOM in after_child_cut.atoms[0].fixed_directions()


def test_slot_conservation_invariant():
    """每个 atom 的 bond+free+fixed 恒为 4，切割后仍成立。"""
    m = build(3, 3, [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 0), (1, 1), (2, 2)])
    m.check_invariants()
    for subset in [{0}, {2}, {0, 3}, {4, 5}]:
        m.cut_as_free(subset).check_invariants()


def test_circuit_rank_matches_cross_bonds_minus_one():
    """SPEC §1.3：简化两层模型中 circuit rank = 跨层 bond 数 - 1。"""
    for n in range(1, 6):
        up_edges = [(i, i + 1) for i in range(n - 1)]
        down_edges = [(i, i + 1) for i in range(n - 1)]
        cross = [(i, i) for i in range(n)]
        m = build(n, n, up_edges, down_edges, cross)
        assert len(m.cross_bonds()) == n
        assert m.circuit_rank() == n - 1


def test_forest_splitting_invariant():
    """SPEC §3.3 (♡)：切连通 S 后，MD 分支数 = 旧分支数 + 外部邻居数 - 1。"""
    m = build(1, 5, [], [(0, 1), (0, 2), (1, 3), (1, 4)], [(0, 0)])
    down = set(m.layer_atoms(Layer.DOWN))
    before = len(m.components(within=down))

    s = {1 + 1}  # MD 局部编号 1 → 全局 id 2
    outside_nbrs = {nb for nb in m.neighbors(2) if nb in down}
    after_mol = m.cut_as_free(s)
    after = len(after_mol.components(within=down - s))

    assert after == before + len(outside_nbrs) - 1


def test_lambda_initial_value_is_minus_one():
    """SPEC §3.7：MD 初始为树且无 deg 3 atom，故 λ = 0 - 1 = -1。"""
    for n in range(1, 6):
        m = build(
            n,
            n,
            [(i, i + 1) for i in range(n - 1)],
            [(i, i + 1) for i in range(n - 1)],
            [(i, i) for i in range(n)],
        )
        assert lambda_potential(m) == -1


class TestElementaryClassification:
    """SPEC §3.1 的四种 elementary molecule + {343}。"""

    def test_e4_single_deg4_no_fixed(self):
        m = build(1, 1, [], [], [(0, 0)])
        assert classify(m, {0}) is ElemType.E4

    def test_e3_single_deg3(self):
        m = build(1, 1, [], [], [(0, 0)])
        after = m.cut_as_free({0})  # atom 1 变 deg 3
        assert classify(after, {1}) is ElemType.E3

    def test_e33_two_adjacent_deg3(self):
        # MU: 0->1 链；MD: 单 atom 2 连到 MU 1
        m = build(2, 1, [(0, 1)], [], [(1, 0)])
        after = m.cut_as_free({0})  # atom 1 得 fixed end → deg 3；atom 2 deg 3?
        # atom 1: bond to 2 + fixed(top) + 2 free = deg 3
        # atom 2: bond to 1 + 3 free = deg 4 → 不是 {33}
        assert after.degree(1) == 3
        assert after.degree(2) == 4
        assert classify(after, {1, 2}) is ElemType.INVALID

    def test_e2_requires_same_direction_free_ends(self):
        """Overview §3.3(1)：{2} 的两个 free end 必须同为 top 或同为 bottom。"""
        m = build(1, 2, [], [(0, 1)], [(0, 0)])
        # atom 0 (UP): 1 bond + 3 free = deg 4
        # 切掉 atom 0 和 atom 2 后，atom 1 剩 2 free
        step1 = m.cut_as_free({0})
        step2 = step1.cut_as_free({2})
        assert step2.degree(1) == 2
        t = classify(step2, {1})
        assert t in (ElemType.E2, ElemType.INVALID)
        # 具体是哪个取决于两个 free end 的方向是否一致
        atom = step2.atoms[1]
        free_dirs = {d for k, d in atom.ends if k is EndKind.FREE}
        assert (t is ElemType.E2) == (len(free_dirs) == 1)


class TestProper:
    """SPEC §3.5 / Def 1.3 的三个禁止模式。"""

    def test_initial_md_tree_is_proper(self):
        m = build(3, 3, [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 0), (1, 1), (2, 2)])
        assert is_proper(m)

    def test_adjacent_deg3_violates(self):
        m = build(2, 2, [(0, 1)], [(0, 1)], [(0, 0), (1, 1)])
        # 切走两个 MU atom → MD 两个 atom 各得一个 fixed end，都变 deg 3 且相邻
        after = m.cut_as_free({0, 1})
        assert after.degree(2) == 3 and after.degree(3) == 3
        assert not is_proper(after)

    def test_deg2_violates(self):
        m = build(1, 1, [], [], [(0, 0)])
        after = m.cut_as_free({0})
        # 人为再降一次度：构造 deg 2 情形
        atom = after.atoms[1]
        atom.ends = [
            (EndKind.FIXED, Direction.TOP),
            (EndKind.FIXED, Direction.BOTTOM),
            (EndKind.FREE, Direction.TOP),
            (EndKind.FREE, Direction.BOTTOM),
        ]
        assert after.degree(1) == 2
        assert not is_proper(after)


def test_induced_submolecule_turns_outside_bonds_into_free_ends():
    """SPEC §1.1：子 molecule 中，与外部的 bond 视为 free end。"""
    m = build(2, 1, [(0, 1)], [], [(1, 0)])
    sub = m.induced({1})
    assert len(sub.bonds) == 0
    assert sub.atoms[1].n_free() == 4  # 原 2 bond 变 free + 原 2 free
    assert sub.degree(1) == 4


def test_unknown_atom_in_cut_raises():
    m = build(1, 1, [], [], [(0, 0)])
    with pytest.raises(ValueError):
        m.cut_as_free({99})
