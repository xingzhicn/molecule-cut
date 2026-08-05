"""transpose_family 的验证（FINDINGS F9 的上界构造）。

不只验证最终数值 #{33} = k，而是**逐条验证证明中的每个断言**——初始度数、
阶段划分、级联触发条件、每个 {343} 打断的跨层 bond 数。这样若证明中某一步
的推理有误，测试会指出具体是哪一步，而不只是「结果不对」。
"""

import pytest

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.elementary import ElemType, classify, is_proper
from molecule_cut.exhaustive import a_of_molecule
from molecule_cut.families import (
    family_for,
    transpose_cross_target,
    transpose_family,
)
from molecule_cut.molecule import Layer

K_RANGE = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12]


def md_local(mol, aid: int) -> int:
    """全局 id → MD 局部编号。"""
    return aid - len(mol.layer_atoms(Layer.UP))


class TestConstruction:
    @pytest.mark.parametrize("k", K_RANGE)
    def test_is_valid_toy_model_one(self, k):
        """Def 1.2：每个 MU atom 恰有一条跨层 bond。"""
        mol = transpose_family(k)
        mol.check_invariants()
        down = set(mol.layer_atoms(Layer.DOWN))
        for u in mol.layer_atoms(Layer.UP):
            assert sum(1 for c in mol.children(u) if c in down) == 1

    @pytest.mark.parametrize("k", K_RANGE)
    def test_cross_map_is_bijection_on_first_3k(self, k):
        """转置映射在 [0,3k) 上是双射；u_{3k} 额外落回 MD_0。"""
        targets = [transpose_cross_target(i, k) for i in range(3 * k)]
        assert sorted(targets) == list(range(3 * k))
        assert transpose_cross_target(3 * k, k) == 0

    @pytest.mark.parametrize("k", K_RANGE)
    def test_block_receives_one_bond_from_each_third(self, k):
        """证明的核心构造性质：MD 块 {3r,3r+1,3r+2} 收到来自
        u_r、u_{k+r}、u_{2k+r} 的三条跨层 bond。"""
        for r in range(k):
            assert transpose_cross_target(r, k) == 3 * r
            assert transpose_cross_target(k + r, k) == 3 * r + 1
            assert transpose_cross_target(2 * k + r, k) == 3 * r + 2

    @pytest.mark.parametrize("k", K_RANGE)
    def test_initial_degrees_all_four(self, k):
        """Def 1.1(3)：初始无 fixed end，全部 deg 4。"""
        mol = transpose_family(k)
        assert all(mol.degree(a) == 4 for a in mol.atoms)
        assert is_proper(mol)


class TestExecution:
    @pytest.mark.parametrize("k", K_RANGE)
    def test_operation_profile_is_rigid(self, k):
        """剖面恒为 (A,B,C,D,cutU) = (k+1, 0, 0, k, 2k)。"""
        rec = run_toy1_algorithm(transpose_family(k))
        assert rec.failed is None, rec.failed
        c = rec.op_counts()
        assert (
            c.get(OpKind.A, 0),
            c.get(OpKind.B, 0),
            c.get(OpKind.CUT33, 0),
            c.get(OpKind.CUT343, 0),
            c.get(OpKind.CUT_U_ONLY, 0),
        ) == (k + 1, 0, 0, k, 2 * k)

    @pytest.mark.parametrize("k", K_RANGE)
    def test_yields_exactly_k_good_components(self, k):
        rec = run_toy1_algorithm(transpose_family(k))
        assert rec.failed is None
        assert rec.n33 == k
        assert rec.n4 == 1

    @pytest.mark.parametrize("k", K_RANGE)
    def test_three_phase_structure(self, k):
        """阶段 1 全是 (a)，阶段 2 全是 {343}，阶段 3 全是 CUT_U。"""
        rec = run_toy1_algorithm(transpose_family(k))
        ops = [s.op for s in rec.steps]
        assert ops[: k + 1] == [OpKind.A] * (k + 1)
        assert ops[k + 1 : 2 * k + 1] == [OpKind.CUT343] * k
        assert ops[2 * k + 1 :] == [OpKind.CUT_U_ONLY] * (2 * k)

    @pytest.mark.parametrize("k", K_RANGE)
    def test_phase1_cuts_top_third_of_mu(self, k):
        """阶段 1 依次切 u_{3k}, u_{3k-1}, …, u_{2k}。"""
        rec = run_toy1_algorithm(transpose_family(k))
        cut = [min(s.subset) for s in rec.steps[: k + 1]]
        assert cut == list(range(3 * k, 2 * k - 1, -1))

    @pytest.mark.parametrize("k", K_RANGE)
    def test_phase2_cuts_consecutive_md_blocks(self, k):
        """阶段 2 的第 r 步切 MD 三元块 {3r, 3r+1, 3r+2}，且确为 {343}。"""
        mol = transpose_family(k)
        rec = run_toy1_algorithm(mol)
        state = mol.copy()
        for s in rec.steps[: k + 1]:
            state = state.cut_as_free(set(s.subset))
        for r, step in enumerate(rec.steps[k + 1 : 2 * k + 1]):
            block = sorted(md_local(mol, a) for a in step.subset)
            assert block == [3 * r, 3 * r + 1, 3 * r + 2]
            assert classify(state, set(step.subset)) is ElemType.E343
            state = state.cut_as_free(set(step.subset))

    @pytest.mark.parametrize("k", K_RANGE)
    def test_each_343_breaks_exactly_two_cross_bonds(self, k):
        """证明的关键计数：块 r 被切时 u_{2k+r} 已在阶段 1 切掉，
        u_r 与 u_{k+r} 仍在，故恰好打断 2 条。"""
        rec = run_toy1_algorithm(transpose_family(k))
        for step in rec.steps[k + 1 : 2 * k + 1]:
            assert step.cross_bonds_broken == 2
        for step in rec.steps[: k + 1]:
            assert step.cross_bonds_broken == 1


class TestIdentitiesAndBounds:
    @pytest.mark.parametrize("k", K_RANGE)
    def test_cutU_identity(self, k):
        """n33 = n − 1 − cutU（由 (1.3)、(1.4) 与 A+B+cutU = n 消元所得）。"""
        n = 3 * k + 1
        rec = run_toy1_algorithm(transpose_family(k))
        cut_u = rec.op_counts().get(OpKind.CUT_U_ONLY, 0)
        assert rec.n33 == n - 1 - cut_u

    @pytest.mark.parametrize("k", K_RANGE)
    def test_satisfies_prop16(self, k):
        """族仍满足 Prop 1.6——它给的是上界，不与已证下界冲突。"""
        n = 3 * k + 1
        rec = run_toy1_algorithm(transpose_family(k))
        assert rec.n33 >= (n - 1) / 5

    @pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
    def test_min_over_all_tiebreaks_is_also_k(self, k):
        """a_n 是 min over tie-breaks；确认没有更低的分支。"""
        assert a_of_molecule(transpose_family(k)).n33_min == k


class TestAllResidueClasses:
    """三个变体合起来覆盖所有 n ≥ 4，给出 a_n ≤ ⌈(n−1)/3⌉。"""

    @pytest.mark.parametrize("n", list(range(4, 40)))
    def test_family_for_hits_the_bound(self, n):
        import math

        mol = family_for(n)
        mol.check_invariants()
        assert len(mol.layer_atoms(Layer.UP)) == n
        rec = run_toy1_algorithm(mol)
        assert rec.failed is None, rec.failed
        assert rec.n33 == math.ceil((n - 1) / 3)
        assert rec.n4 == 1

    @pytest.mark.parametrize("n", list(range(4, 40)))
    def test_family_for_is_valid_toy1(self, n):
        mol = family_for(n)
        down = set(mol.layer_atoms(Layer.DOWN))
        for u in mol.layer_atoms(Layer.UP):
            assert sum(1 for c in mol.children(u) if c in down) == 1

    @pytest.mark.parametrize("n", [5, 6, 8, 9, 11, 12])
    def test_min_over_tiebreaks_matches(self, n):
        import math

        assert a_of_molecule(family_for(n)).n33_min == math.ceil((n - 1) / 3)


class TestVariantExecutionMatchesProof:
    """两个变体的操作序列逐步匹配 FINDINGS F9 中写出的证明。

    只断言最终 #{33} 不足以验证证明——证明声称了**具体的执行序列**，
    这里把该序列本身作为断言，任一阶段错位都会被指出。
    """

    @pytest.mark.parametrize("k", [1, 2, 3, 4, 5, 6, 7, 8])
    def test_padded_execution(self, k):
        """(k+1) 次 (a) → k 次 {343} → 2k 次 CUT_U → 1 次 (b)。"""
        from molecule_cut.families import padded_family

        rec = run_toy1_algorithm(padded_family(k))
        assert rec.failed is None, rec.failed
        assert [s.op for s in rec.steps] == (
            [OpKind.A] * (k + 1) + [OpKind.CUT343] * k + [OpKind.CUT_U_ONLY] * (2 * k) + [OpKind.B]
        )
        assert rec.n33 == k + 1

    @pytest.mark.parametrize("k", [2, 3, 4, 5, 6, 7, 8, 10, 12])
    def test_trimmed_execution(self, k):
        """k 次 (a) → (k−1) 次 {343} → 1 次 (a) → 1 次 {33} → (2k−1) 次 CUT_U。"""
        from molecule_cut.families import trimmed_family

        rec = run_toy1_algorithm(trimmed_family(k))
        assert rec.failed is None, rec.failed
        assert [s.op for s in rec.steps] == (
            [OpKind.A] * k
            + [OpKind.CUT343] * (k - 1)
            + [OpKind.A, OpKind.CUT33]
            + [OpKind.CUT_U_ONLY] * (2 * k - 1)
        )
        assert rec.n33 == k

    @pytest.mark.parametrize("k", [2, 3, 4, 5, 6, 8, 12])
    def test_trimmed_final_cut33_breaks_only_one(self, k):
        """证明第 5 步的关键：尾对的 {33} 只断 1 条，故 Σb_i = 2g − 1。"""
        from molecule_cut.families import trimmed_family

        rec = run_toy1_algorithm(trimmed_family(k))
        breaks = [s.cross_bonds_broken for s in rec.steps if s.op in (OpKind.CUT343, OpKind.CUT33)]
        assert breaks == [2] * (k - 1) + [1]
        assert sum(breaks) == 2 * len(breaks) - 1

    @pytest.mark.parametrize("k", [1, 2, 3, 4, 5, 6, 8])
    def test_padded_all_cleanups_break_two(self, k):
        """padded 的 k 次 {343} 各断 2 条，故 Σb_i = 2g 取等号。"""
        from molecule_cut.families import padded_family

        rec = run_toy1_algorithm(padded_family(k))
        breaks = [s.cross_bonds_broken for s in rec.steps if s.op in (OpKind.CUT343, OpKind.CUT33)]
        assert breaks == [2] * k
