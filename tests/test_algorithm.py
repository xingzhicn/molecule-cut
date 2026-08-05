"""Def 1.4 算法与 Prop 1.6 的复现测试（SPEC §3.6, §3.7）。

最重要的是 `TestProp16Identities`：它把讲义证明 part 4 的三条计数关系
写成断言。若实现与讲义有任何语义偏差，这些恒等式几乎必然被打破。
"""

import pytest

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.builders import build
from molecule_cut.molecule import Layer

# SPEC §3.7：每种操作最多能打断的跨层 bond 数
MAX_CROSS_BROKEN = {
    OpKind.A: 1,
    OpKind.B: 2,
    OpKind.CUT33: 3,
    OpKind.CUT343: 4,
}

# SPEC §3.7：每种操作对 λ 的改变量
LAMBDA_DELTA = {
    OpKind.A: 1,
    OpKind.B: 0,
    OpKind.CUT33: -1,
    OpKind.CUT343: -1,
    OpKind.CUT3: 0,
}


def chain_toy1(n: int, md_edges=None):
    """|MU| = n 的链式 toy model I：MU 为链（唯一最低 atom），MD 默认也是链。"""
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = md_edges if md_edges is not None else [(i, i + 1) for i in range(n - 1)]
    cross = [(i, i) for i in range(n)]
    return build(n, n, up_edges, down_edges, cross)


def star_md_toy1(n: int):
    """MD 为星形（0 为根，其余为其 child）。注意每 atom 至多 2 child，
    故 n > 3 时退化为二叉树形。"""
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = []
    for i in range(1, n):
        parent = (i - 1) // 2
        down_edges.append((parent, i))
    cross = [(i, i) for i in range(n)]
    return build(n, n, up_edges, down_edges, cross)


class TestProp16Identities:
    """Prop 1.6 证明 part 3-4 的计数关系。"""

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_n4_is_exactly_one(self, n, builder):
        """Prop 1.6：#{4} = 1（MU 恰有一个最低 atom，且只有它可能是 deg 4）。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        assert rec.n4 == 1

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_n33_lower_bound(self, n, builder):
        """Prop 1.6：#{33} ≥ (|MU| - 1) / 5。这是 M1 的正确性 Gate。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        assert rec.n33 >= (n - 1) / 5

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_identity_1_3(self, n, builder):
        """(1.3)：#{33} = #(b) + #cut{33} + #cut{343}。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        c = rec.op_counts()
        rhs = c.get(OpKind.B, 0) + c.get(OpKind.CUT33, 0) + c.get(OpKind.CUT343, 0)
        assert rec.n33 == rhs

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_identity_1_4(self, n, builder):
        """(1.4)：#(a) - (#cut{33} + #cut{343}) = 1，源自 λ 从 -1 走到 0。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        c = rec.op_counts()
        lhs = c.get(OpKind.A, 0) - (c.get(OpKind.CUT33, 0) + c.get(OpKind.CUT343, 0))
        assert lhs == 1

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_bond_counting_inequality(self, n, builder):
        """证明 part 4 的第三条约束：
        #(a) + 2#(b) + 3#cut{33} + 4#cut{343} ≥ |MU|。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        c = rec.op_counts()
        lhs = (
            c.get(OpKind.A, 0)
            + 2 * c.get(OpKind.B, 0)
            + 3 * c.get(OpKind.CUT33, 0)
            + 4 * c.get(OpKind.CUT343, 0)
        )
        assert lhs >= n


class TestStepwiseInvariants:
    """逐步不变量：λ 变化量与跨层 bond 打断数的上界（SPEC §3.7 表）。"""

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_lambda_delta_per_operation(self, n, builder):
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        for s in rec.steps:
            if s.op in LAMBDA_DELTA:
                assert s.lam_after - s.lam_before == LAMBDA_DELTA[s.op], (
                    f"op {s.op.value} changed lambda by "
                    f"{s.lam_after - s.lam_before}, expected {LAMBDA_DELTA[s.op]}"
                )

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_cross_bond_breaking_bounds(self, n, builder):
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        for s in rec.steps:
            if s.op in MAX_CROSS_BROKEN:
                assert s.cross_bonds_broken <= MAX_CROSS_BROKEN[s.op], (
                    f"op {s.op.value} broke {s.cross_bonds_broken} cross bonds, "
                    f"max allowed {MAX_CROSS_BROKEN[s.op]}"
                )

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_lambda_starts_at_minus_one_ends_at_zero(self, n, builder):
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        assert rec.steps[0].lam_before == -1
        assert rec.steps[-1].lam_after == 0

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_all_atoms_consumed(self, n, builder):
        """算法结束时每个 atom 恰好被切一次。"""
        mol = builder(n)
        rec = run_toy1_algorithm(mol)
        assert rec.failed is None, rec.failed
        cut = [a for s in rec.steps for a in s.subset]
        assert sorted(cut) == sorted(mol.atoms)

    @pytest.mark.parametrize("n", range(1, 8))
    @pytest.mark.parametrize("builder", [chain_toy1, star_md_toy1])
    def test_every_component_is_elementary(self, n, builder):
        """Prop 1.6 part 1-2：过程中只产出 elementary 组件。"""
        rec = run_toy1_algorithm(builder(n))
        assert rec.failed is None, rec.failed
        # 每个 {33} 覆盖 2 个 atom，其余各覆盖 1 个。{343} 已被拆记为
        # 一个 {33} + 一个 {3}（共 3 个 atom），故无需额外修正。
        total = rec.n33 * 2 + rec.n4 + rec.n3 + rec.n2
        assert total == 2 * n, (
            f"atom count mismatch: components cover {total}, molecule has {2 * n}"
        )


class TestScore:
    @pytest.mark.parametrize("n", range(1, 8))
    def test_score_formula(self, n):
        rec = run_toy1_algorithm(chain_toy1(n))
        assert rec.failed is None
        assert rec.score(d=3) == rec.n33 - 30 * rec.n4

    def test_toy1_structural_precondition(self):
        """Toy I 要求 MU 恰有一个最低 atom（Prop 1.6 part 3）。"""
        m = chain_toy1(4)
        up = set(m.layer_atoms(Layer.UP))
        lowest = [a for a in up if not (set(m.children(a)) & up)]
        assert len(lowest) == 1
