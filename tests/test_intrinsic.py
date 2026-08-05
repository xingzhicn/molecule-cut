"""内在量 v_d(M) 与结构性观察的测试。

Observation A（可证，讲义未显式指出）：
  在任何简化两层模型中（Def 1.1(3)：初始无 fixed end），所有 atom 初始
  deg 4。而 {33}/{3}/{2} 均要求至少一个 fixed end，故**任何合法切割序列的
  第一步必须切出一个 {4}**。因此 #{4} ≥ 1 对所有算法成立，Prop 1.6 得到的
  #{4} = 1 在这一项上是最优的。
"""

import math

import pytest

from molecule_cut.algorithm import run_toy1_algorithm
from molecule_cut.builders import build
from molecule_cut.elementary import ElemType, classify
from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.exhaustive import _legal_moves, v_d, w, w_reference


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_observation_a_first_cut_must_be_e4(n):
    """初始状态下唯一可切的 elementary 组件类型是 {4}。"""
    for mol in enumerate_toy1(n):
        moves = _legal_moves(mol)
        assert moves, f"no legal move on initial molecule (n={n})"
        types = {classify(mol, set(mv)) for mv in moves}
        assert types == {ElemType.E4}, (
            f"n={n}: initial legal cuts should be exactly {{4}}, got {types}"
        )


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_all_atoms_initially_degree_four(n):
    """Def 1.1(3) 的直接推论：初始无 fixed end ⟹ 全部 deg 4。"""
    for mol in enumerate_toy1(n):
        for aid in mol.atoms:
            assert mol.atoms[aid].n_fixed() == 0
            assert mol.degree(aid) == 4


def test_w_is_not_unconditionally_v_d_plus_10d():
    """W 最大化好分量，v_d 最大化带 {4} 惩罚的得分，目标不可混同。"""
    mol = build(3, 3, [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(2, 2)])
    assert w_reference(mol) == 1
    assert w(mol) == 1
    assert v_d(mol, d=3) == -30
    assert v_d(mol, d=3) + 30 == 0


def test_v_d_monotone_in_d():
    """d 增大时 v_d 不增（惩罚项变重）。"""
    for mol in enumerate_toy1(3):
        v2, v3, v4 = v_d(mol, d=2), v_d(mol, d=3), v_d(mol, d=4)
        assert v2 >= v3 >= v4


@pytest.mark.parametrize("d", [2, 3])
@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_intrinsic_vd_dominates_sharp_benchmark_witness(n, d):
    """Regression for (T_vd): the intrinsic maximum dominates the witness run."""
    bound = math.ceil((n - 1) / 3) - 10 * d
    for mol in enumerate_toy1(n):
        witness = run_toy1_algorithm(mol)
        assert witness.failed is None, witness.failed
        assert witness.n4 == 1
        assert witness.n33 >= math.ceil((n - 1) / 3)

        intrinsic = v_d(mol, d=d)
        assert intrinsic is not None
        assert intrinsic >= witness.score(d)
        assert intrinsic >= bound
