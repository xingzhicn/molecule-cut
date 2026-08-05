"""讲义计数系数可达性的反例见证（FINDINGS F3 修订）。

背景：在 |MU| ≤ 6 的全枚举中，cut{33} 最多打断 2 条跨层 bond、cut{343} 最多
3 条，低于讲义 Prop 1.6 证明 part 4 使用的系数 3 与 4。这一度让人猜测系数可
收紧、从而改进常数 1/5。

**该猜想是错的。** 在更大的 |MU| 上两个系数都能达到。下面固化具体见证，
防止「小样本假象」被再次当成定理。见证由随机搜索发现，此处只做确定性重放。
"""

import pytest

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.elementary import ElemType, classify
from molecule_cut.molecule import Layer
from molecule_cut.serialize import from_json

# cut{33} 打断 3 条跨层 bond：|MU|=8, |MD|=7，默认 tie-breaking 下第 6 步
WITNESS_CUT33_BREAKS_3 = (
    '{"cross":[[0,1],[1,0],[2,3],[3,0],[4,5],[5,2],[6,4],[7,6]],"down":7,'
    '"down_edges":[[0,1],[0,2],[1,5],[2,3],[2,4],[5,6]],"up":8,'
    '"up_edges":[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7]]}'
)

# cut{343} 打断 4 条跨层 bond：|MU|=11, |MD|=13，默认 tie-breaking 下第 6 步
WITNESS_CUT343_BREAKS_4 = (
    '{"cross":[[0,11],[1,3],[2,0],[3,0],[4,1],[5,9],[6,2],[7,6],[8,7],[9,12],[10,4]],'
    '"down":13,"down_edges":[[0,1],[0,3],[1,2],[1,5],[2,7],[3,4],[3,6],[5,8],[5,9],'
    '[7,10],[7,11],[11,12]],"up":11,'
    '"up_edges":[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]]}'
)


def _is_toy1(mol) -> bool:
    """Def 1.2：每个 MU atom 恰有一条跨层 bond。"""
    down = set(mol.layer_atoms(Layer.DOWN))
    return all(sum(1 for c in mol.children(u) if c in down) == 1 for u in mol.layer_atoms(Layer.UP))


def _replay_until(mol, step_index: int):
    """把 molecule 重放到第 step_index 步之前的状态。"""
    rec = run_toy1_algorithm(mol)
    assert rec.failed is None, rec.failed
    cur = mol.copy()
    for s in rec.steps[:step_index]:
        cur = cur.cut_as_free(set(s.subset))
    return rec, cur


def _cross_bonds_broken(mol, subset: set[int]) -> list[tuple[int, int]]:
    """独立重算：被切走的 subset 会打断哪些当前还存在的跨层 bond。"""
    return [
        (b.parent, b.child)
        for b in mol.cross_bonds()
        if (b.parent in subset) or (b.child in subset)
    ]


@pytest.mark.parametrize(
    "witness,op,expected_broken,step_index,expected_type",
    [
        (WITNESS_CUT33_BREAKS_3, OpKind.CUT33, 3, 6, ElemType.E33),
        (WITNESS_CUT343_BREAKS_4, OpKind.CUT343, 4, 6, ElemType.E343),
    ],
    ids=["cut33_breaks_3", "cut343_breaks_4"],
)
def test_notes_coefficients_are_attained(witness, op, expected_broken, step_index, expected_type):
    """讲义系数可达：独立重算打断数，并核对该步切出的确实是对应类型。"""
    mol = from_json(witness)
    mol.check_invariants()
    assert _is_toy1(mol), "witness must be a toy model I instance"

    rec, state = _replay_until(mol, step_index)
    step = rec.steps[step_index]
    subset = set(step.subset)

    assert step.op is op
    assert classify(state, subset) is expected_type
    for aid in subset:
        assert state.degree(aid) in (3, 4)

    broken = _cross_bonds_broken(state, subset)
    assert len(broken) == expected_broken, (
        f"expected {expected_broken} cross bonds broken, recomputed {broken}"
    )
    assert step.cross_bonds_broken == expected_broken


def test_witnesses_satisfy_prop16():
    """两个见证本身仍满足 Prop 1.6（系数可达不等于定理失效）。"""
    for witness in (WITNESS_CUT33_BREAKS_3, WITNESS_CUT343_BREAKS_4):
        mol = from_json(witness)
        n = len(mol.layer_atoms(Layer.UP))
        rec = run_toy1_algorithm(mol)
        assert rec.failed is None
        assert rec.n4 == 1
        assert rec.n33 >= (n - 1) / 5
