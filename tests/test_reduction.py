"""旧 `P` 路线的回归检查与已验证不变量。

本文件保留对充分条件 `Σ b_i ≤ 2g` 的历史回归测试；该条件不是当前
`P*` 流证明的必要目标。恒等式失败说明实现有 bug；不等式失败仍是一个
有价值的反例记录，但不再推翻已证的 benchmark 定理。
"""

import random

import pytest

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.builders import build
from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.families import family_for
from molecule_cut.molecule import Direction, Layer

CLEANUP_OPS = (OpKind.CUT33, OpKind.CUT343)
SEED = 20260730


def random_toy1(n: int, nd: int, rng: random.Random):
    edges, child_count = [], {0: 0}
    for v in range(1, nd):
        avail = [u for u in range(v) if child_count[u] < 2]
        p = rng.choice(avail)
        edges.append((p, v))
        child_count[p] += 1
        child_count[v] = 0
    cap = [2] * nd
    for _, c in edges:
        cap[c] -= 1
    slots = [j for j in range(nd) for _ in range(cap[j])]
    if len(slots) < n:
        return None
    cross = list(enumerate(rng.sample(slots, n)))
    try:
        return build(n, nd, [(i, i + 1) for i in range(n - 1)], edges, cross)
    except ValueError:
        return None


def cleanup_stats(rec) -> tuple[int, int]:
    """返回 (g, Σb_i)：清理操作个数与它们打断的跨层 bond 总数。"""
    steps = [s for s in rec.steps if s.op in CLEANUP_OPS]
    return len(steps), sum(s.cross_bonds_broken for s in steps)


class TestProvenIdentity:
    """[已证] n33 = n − 1 − #cutU，由 (1.3)、(1.4) 与 A+B+cutU = n 消元。"""

    @pytest.mark.parametrize("n", [2, 3, 4, 5])
    def test_cutU_identity_enumerated(self, n):
        for mol in enumerate_toy1(n):
            rec = run_toy1_algorithm(mol)
            assert rec.failed is None
            cut_u = rec.op_counts().get(OpKind.CUT_U_ONLY, 0)
            assert rec.n33 == n - 1 - cut_u

    @pytest.mark.parametrize("n", list(range(6, 20)))
    def test_cutU_identity_random(self, n):
        rng = random.Random(SEED + n)
        checked = 0
        for _ in range(40):
            mol = random_toy1(n, rng.randint(max(2, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            checked += 1
            cut_u = rec.op_counts().get(OpKind.CUT_U_ONLY, 0)
            assert rec.n33 == n - 1 - cut_u
        assert checked > 0, f"n={n}: 没有生成任何有效实例，测试空转"


class TestVerifiedStructuralLemma:
    """[已验证，未证明] 携带 2 条跨层 bond 的 MD atom 的结构约束。

    观察：若 MD atom a 在被切时携带 2 条跨层 bond，则 a 无 MD 父（两个 parent
    槽都被跨层 bond 占用）；若 a 还是 deg 3，则其 fixed end 必为 bottom，
    即 a 的某个 child 更早被移除。
    """

    @pytest.mark.parametrize("n", list(range(6, 16)))
    def test_two_cross_bond_atoms_have_no_md_parent(self, n):
        rng = random.Random(SEED * 2 + n)
        seen = 0
        for _ in range(40):
            mol = random_toy1(n, rng.randint(max(2, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            for step in rec.steps:
                if step.op in CLEANUP_OPS:
                    up = set(state.layer_atoms(Layer.UP))
                    for a in step.subset:
                        if sum(1 for p in state.parents(a) if p in up) == 2:
                            seen += 1
                            assert not [p for p in state.parents(a) if p not in up]
                            if state.degree(a) == 3:
                                assert Direction.BOTTOM in state.atoms[a].fixed_directions()
                state = state.cut_as_free(set(step.subset))
        # 不断言 seen > 0：小 n 上可能不出现，这是允许的
        assert seen >= 0


class TestUnprovenSharpInequality:
    """[未证明且非必要] Σ b_i ≤ 2g —— 旧路线的充分条件（见 FINDINGS F10）。

    当前 benchmark 下界由带 `B` 修正项的 `P*` 流证明给出；若本测试失败，
    只说明旧的更强条件 `P` 有反例，不要把它当作当前定理的失败。
    """

    @pytest.mark.parametrize("n", [2, 3, 4, 5])
    def test_sharp_inequality_enumerated(self, n):
        for mol in enumerate_toy1(n):
            rec = run_toy1_algorithm(mol)
            g, total = cleanup_stats(rec)
            assert total <= 2 * g, f"n={n}: Σb_i={total} > 2g={2 * g} —— F7 的反例！"

    @pytest.mark.parametrize("n", list(range(6, 20)))
    def test_sharp_inequality_random(self, n):
        rng = random.Random(SEED * 3 + n)
        checked = 0
        for _ in range(60):
            mol = random_toy1(n, rng.randint(max(2, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            checked += 1
            g, total = cleanup_stats(rec)
            assert total <= 2 * g, f"n={n}: Σb_i={total} > 2g={2 * g} —— F7 的反例！"
        assert checked > 0, f"n={n}: 没有生成任何有效实例，测试空转"

    @pytest.mark.parametrize("n", list(range(4, 30)))
    def test_families_are_near_tight(self, n):
        """F9 的族把不等式压到紧或差 1，说明 2 这个系数是尖锐的。

        n ≡ 1, 2 (mod 3) 时取等号 Σb_i = 2g；n ≡ 0 (mod 3) 时
        `trimmed_family` 给出 Σb_i = 2g − 1，差额被 ⌈·⌉ 吸收，
        故仍达到目标 n33 = ⌈(n−1)/3⌉。
        """
        rec = run_toy1_algorithm(family_for(n))
        g, total = cleanup_stats(rec)
        assert g > 0
        assert total <= 2 * g
        assert total >= 2 * g - 1, f"n={n}: Σb_i={total} 距 2g={2 * g} 差 {2 * g - total}"
        if n % 3 != 0:
            assert total == 2 * g

    def test_inequality_is_sufficient_legacy_route(self):
        """验证旧归约链：Σb_i ≤ 2g ⟹ n33 ≥ (n−1)/3。"""
        import math

        rng = random.Random(SEED * 5)
        for n in range(6, 18):
            for _ in range(25):
                mol = random_toy1(n, rng.randint(max(2, n - 2), n + 3), rng)
                if mol is None:
                    continue
                rec = run_toy1_algorithm(mol, check=False)
                if rec.failed:
                    continue
                g, total = cleanup_stats(rec)
                if total <= 2 * g:
                    assert rec.n33 >= math.ceil((n - 1) / 3) or rec.n33 >= (n - 1) / 3


class TestStructuralLemmasAreFalseInGeneral:
    """[反例] 两条「结构引理」作为任意 molecule 的命题是**假的**。

    它们在算法实际运行中从未被违反（见 FINDINGS F10），但这依赖**可达性**，
    不是局部结构性质。下面的显式反例把这一区别钉死，防止再次误以为
    它们可由槽位计数证明。

    共同反例：MD 中一个 atom 同时有两个 MD 父。此时「每个分支恰有一个根」
    的论证失效——MD 的定向不是有根树，一个 atom 可以有 2 个 parent 槽被
    MD 边占用，因而可以有多个 atom 都没有 MD 父。
    """

    @staticmethod
    def two_parent_molecule():
        """MD: atom a → b ← c，即 b 有两个 MD 父；a、c 各挂 2 条跨层 bond。"""
        return build(
            4, 3, [(0, 1), (1, 2), (2, 3)], [(0, 1), (2, 1)], [(0, 0), (1, 0), (2, 2), (3, 2)]
        )

    def test_two_atoms_can_carry_two_cross_bonds_each(self):
        """反驳「至多一个 atom 携带 2 条跨层 bond」的结构版本。"""
        mol = self.two_parent_molecule()
        mol.check_invariants()
        up = set(mol.layer_atoms(Layer.UP))
        carriers = [
            a for a in mol.layer_atoms(Layer.DOWN) if sum(1 for p in mol.parents(a) if p in up) == 2
        ]
        assert len(carriers) == 2, "应存在两个各携带 2 条跨层 bond 的 MD atom"
        for a in carriers:
            assert not [p for p in mol.parents(a) if p not in up], "二者都应无 MD 父"

    def test_one_removal_can_create_two_bottom_fixed_ends(self):
        """反驳「一次移除至多产生一个 bottom fixed end」的结构版本。"""
        mol = self.two_parent_molecule()
        down = set(mol.layer_atoms(Layer.DOWN))
        target = {5}  # 有两个 MD 父的那个 atom
        before = {a: mol.atoms[a].fixed_directions() for a in down - target}
        after = mol.cut_as_free(target)
        new_bottom = [
            a
            for a in down - target
            if Direction.BOTTOM in after.atoms[a].fixed_directions()
            and Direction.BOTTOM not in before[a]
        ]
        assert len(new_bottom) == 2, f"应新增 2 个 bottom fixed end，实得 {new_bottom}"

    def test_lemmas_still_hold_along_actual_runs(self):
        """但沿算法实际执行，两条观察仍未被违反——差别在可达性。"""
        rng = random.Random(SEED * 7)
        checked_a = checked_b = 0
        for n in range(6, 14):
            for _ in range(15):
                mol = random_toy1(n, rng.randint(max(2, n - 2), n + 3), rng)
                if mol is None:
                    continue
                rec = run_toy1_algorithm(mol, check=False)
                if rec.failed:
                    continue
                state = mol.copy()
                for step in rec.steps:
                    sub = set(step.subset)
                    up = set(state.layer_atoms(Layer.UP))
                    down = set(state.layer_atoms(Layer.DOWN))
                    if step.op in CLEANUP_OPS:
                        checked_a += 1
                        carriers = sum(
                            1 for a in sub if sum(1 for p in state.parents(a) if p in up) == 2
                        )
                        assert carriers <= 1
                    if sub & down:
                        checked_b += 1
                        before = {a: state.atoms[a].fixed_directions() for a in down - sub}
                        after = state.cut_as_free(sub)
                        new_bottom = sum(
                            1
                            for a in down - sub
                            if Direction.BOTTOM in after.atoms[a].fixed_directions()
                            and Direction.BOTTOM not in before[a]
                        )
                        assert new_bottom <= 1
                    state = state.cut_as_free(sub)
        assert checked_a > 0 and checked_b > 0, "测试空转"


def random_toy1_full(n: int, nd: int, rng: random.Random):
    """MD 定向扩展子域：允许 atom 有两个 MD 父，但 MU 仍固定为链。

    这能测试早期的 MD 单父盲区，却不是完整 Toy I 定义域；完整小规模审计由
    `enumerate_toy1_full_labeled` 负责，见 FINDINGS F11。
    """
    if nd == 1:
        edges = []
    elif nd == 2:
        edges = [(0, 1)]
    else:
        seq = [rng.randrange(nd) for _ in range(nd - 2)]
        deg = [1] * nd
        for x in seq:
            deg[x] += 1
        edges = []
        for x in seq:
            for leaf in range(nd):
                if deg[leaf] == 1:
                    edges.append((leaf, x))
                    deg[leaf] -= 1
                    deg[x] -= 1
                    break
        rest = [v for v in range(nd) if deg[v] == 1]
        edges.append((rest[0], rest[1]))

    oriented, npar, nchi = [], [0] * nd, [0] * nd
    for u, v in edges:
        if rng.random() < 0.5:
            u, v = v, u
        if nchi[u] >= 2 or npar[v] >= 2:
            u, v = v, u
            if nchi[u] >= 2 or npar[v] >= 2:
                return None
        oriented.append((u, v))
        nchi[u] += 1
        npar[v] += 1

    slots = [j for j in range(nd) for _ in range(2 - npar[j])]
    if len(slots) < n:
        return None
    try:
        return build(
            n,
            nd,
            [(i, i + 1) for i in range(n - 1)],
            oriented,
            list(enumerate(rng.sample(slots, n))),
        )
    except ValueError:
        return None


class TestSharpInequalityOnCorrectedSpace:
    """命题 P 在 MD 定向扩展子域上的检验（FINDINGS F11 修正后）。

    此前 `TestUnprovenSharpInequality` 用的 `random_toy1` 与枚举器共享盲区，
    只覆盖 MD 为单父有根树的子类。这里改用允许双父、但仍固定 MU 链的生成器重测。
    """

    @pytest.mark.parametrize("n", list(range(6, 18)))
    def test_sum_b_at_most_2g(self, n):
        rng = random.Random(90210 + n)
        checked = two_parent = 0
        for _ in range(80):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            checked += 1
            down = set(mol.layer_atoms(Layer.DOWN))
            if any(sum(1 for p in mol.parents(a) if p in down) == 2 for a in down):
                two_parent += 1
            g, total = cleanup_stats(rec)
            assert total <= 2 * g, f"n={n}: Σb_i={total} > 2g={2 * g} —— F7 的反例！"
        assert checked > 0, f"n={n}: 测试空转"
        assert two_parent > 0, f"n={n}: 修正后的生成器仍未产出双父实例，盲区未修好"


class TestBothLemmasFailOnCorrectedSpace:
    """记录：两条引理在完整定义域上**沿算法运行也不成立**（F11）。

    F10 曾据受限空间的零违反统计称二者「是可达性性质」；该说法已作废。
    """

    def test_lemma_a_violated_along_runs(self):
        """存在清理操作，其切割集内两个 atom 各带 2 条跨层 bond。"""
        rng = random.Random(555)
        for _ in range(20000):
            n = rng.randint(6, 16)
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            for step in rec.steps:
                if step.op in CLEANUP_OPS:
                    up = set(state.layer_atoms(Layer.UP))
                    carriers = sum(
                        1 for a in step.subset if sum(1 for p in state.parents(a) if p in up) == 2
                    )
                    if carriers >= 2:
                        return  # 找到违反实例，引理 A 作为可达性命题为假
                state = state.cut_as_free(set(step.subset))
        pytest.fail("未能在修正空间找到引理 A 的违反实例；F11 的记录需复核")

    def test_lemma_b_violated_along_runs(self):
        """存在一次 MD 移除产生多于一个 bottom fixed end。"""
        rng = random.Random(31415)
        for _ in range(20000):
            n = rng.randint(6, 16)
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            for step in rec.steps:
                sub = set(step.subset)
                down = set(state.layer_atoms(Layer.DOWN))
                if sub & down:
                    before = {a: state.atoms[a].fixed_directions() for a in down - sub}
                    after = state.cut_as_free(sub)
                    new_bottom = sum(
                        1
                        for a in down - sub
                        if Direction.BOTTOM in after.atoms[a].fixed_directions()
                        and Direction.BOTTOM not in before[a]
                    )
                    if new_bottom >= 2:
                        return
                state = state.cut_as_free(sub)
        pytest.fail("未能在修正空间找到引理 B 的违反实例；F11 的记录需复核")


class TestProvedLemmaC:
    """[已证] 携带 2 条跨层 bond 的 deg-3 MD atom，其 fixed end 必在 child 槽。

    证明（纯槽位计数，**不假设 MD 是有根树**——这正是引理 A/B 出错之处）：
    跨层 bond 占 parent 槽（Def 1.1(2)：MU atom 是 MD atom 的父）。该 atom 的
    两个 parent 槽都被跨层 bond 占满；deg 3 表示 bond + free = 3，而总槽位为 4，
    故恰有 1 个 fixed end；它不可能落在 parent 槽，只能落在 child 槽。∎

    推论：这样的 atom 必定更早失去过一个 child——即存在一次更早的移除事件。
    这是摊还论证的付费方。
    """

    @pytest.mark.parametrize("n", list(range(6, 16)))
    def test_two_cross_bond_deg3_has_bottom_fixed_end(self, n):
        rng = random.Random(6060 + n)
        seen = 0
        for _ in range(60):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            for step in rec.steps:
                up = set(state.layer_atoms(Layer.UP))
                for a in set(state.layer_atoms(Layer.DOWN)):
                    if state.degree(a) != 3:
                        continue
                    if sum(1 for p in state.parents(a) if p in up) == 2:
                        seen += 1
                        assert Direction.BOTTOM in state.atoms[a].fixed_directions()
                        assert Direction.TOP not in state.atoms[a].fixed_directions()
                state = state.cut_as_free(set(step.subset))
        assert seen >= 0


class TestPrefixProperties:
    """[未证] 两条比命题 P 更强的**前缀**性质——摊还论证的目标形式。

    P1  对每个 j：Σ_{i≤j} b_i ≤ 2j     （命题 P 是其 j = g 的特例）
    P2  对每个 j：A_j ≥ j + 1          （恒等式 (1.4) 的局部版本，等价于 λ ≥ 1）

    P1 成立意味着存在一个逐步维持的势 Φ_j = 2j − Σ_{i≤j} b_i ≥ 0，
    比全局求和更适合证明。
    """

    @pytest.mark.parametrize("n", list(range(6, 17)))
    def test_prefix_sum_bounded(self, n):
        rng = random.Random(777 + n)
        checked = 0
        for _ in range(60):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            checked += 1
            running = 0
            for j, step in enumerate((s for s in rec.steps if s.op in CLEANUP_OPS), start=1):
                running += step.cross_bonds_broken
                assert running <= 2 * j, f"n={n}: 前缀 j={j} 处 Σb={running} > {2 * j}"
        assert checked > 0, f"n={n}: 测试空转"

    @pytest.mark.parametrize("n", list(range(6, 17)))
    def test_a_ops_lead_cleanups(self, n):
        rng = random.Random(2024 + n)
        checked = 0
        for _ in range(60):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            checked += 1
            a_count = cleanup_count = 0
            for step in rec.steps:
                if step.op is OpKind.A:
                    a_count += 1
                elif step.op in CLEANUP_OPS:
                    cleanup_count += 1
                    assert a_count >= cleanup_count, (
                        f"n={n}: 第 {cleanup_count} 次清理时只有 {a_count} 次 (a)"
                    )
        assert checked > 0, f"n={n}: 测试空转"


class TestProvedLemmaD:
    """[已证] 逐项界：{33} 的 b_i ≤ 3，{343} 的 b_i ≤ 4。

    {33}：两原子由 bond 相连，作为子的那个有一个 parent 槽被占，故至多带 1 条
    跨层 bond；另一个至多 2 条。合计 ≤ 3。
    {343}：按三种定向逐一计数，三种都得 ≤ 4（且都能取到）。

    证明只用槽位计数与相邻关系，不假设 MD 是有根树。
    """

    @pytest.mark.parametrize("n", list(range(6, 16)))
    def test_per_operation_bounds(self, n):
        rng = random.Random(1357 + n)
        seen33 = seen343 = 0
        for _ in range(50):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            for step in rec.steps:
                if step.op is OpKind.CUT33:
                    seen33 += 1
                    assert step.cross_bonds_broken <= 3
                elif step.op is OpKind.CUT343:
                    seen343 += 1
                    assert step.cross_bonds_broken <= 4
        assert seen33 + seen343 > 0, f"n={n}: 测试空转"

    @pytest.mark.parametrize("n", list(range(6, 16)))
    def test_cut33_child_carries_at_most_one(self, n):
        """{33} 界的证明依据：作为子的那个原子至多带 1 条跨层 bond。"""
        rng = random.Random(2468 + n)
        for _ in range(40):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            for step in rec.steps:
                if step.op is OpKind.CUT33:
                    up = set(state.layer_atoms(Layer.UP))
                    a, b = sorted(step.subset)
                    child = b if a in state.parents(b) else a
                    assert sum(1 for p in state.parents(child) if p in up) <= 1
                state = state.cut_as_free(set(step.subset))


class TestChargingSchemeFails:
    """[负结果] 记录一个已尝试并失败的 charging 方案，避免重复投入。

    方案：由引理 C，每个双跨层 bond 的 deg-3 原子都有 bottom fixed end，
    把它造成的透支记到「造成该 bottom fixed end 的事件」上。

    失败：该事件可能是 b = 2 的清理（无结余），也可能是操作 (b)（不在 Φ 内）。
    """

    def test_paying_event_can_lack_slack(self):
        rng = random.Random(24680)
        found_no_slack = False
        for _ in range(4000):
            n = rng.randint(6, 15)
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            rec = run_toy1_algorithm(mol, check=False)
            if rec.failed:
                continue
            state = mol.copy()
            bottom_source: dict[int, int] = {}
            for idx, step in enumerate(rec.steps):
                sub = set(step.subset)
                up = set(state.layer_atoms(Layer.UP))
                down = set(state.layer_atoms(Layer.DOWN))
                if step.op in CLEANUP_OPS:
                    for a in sub:
                        if (
                            state.degree(a) == 3
                            and sum(1 for p in state.parents(a) if p in up) == 2
                        ):
                            src = bottom_source.get(a)
                            if src is not None:
                                ev = rec.steps[src]
                                if ev.op not in CLEANUP_OPS or ev.cross_bonds_broken >= 2:
                                    found_no_slack = True
                before = {a: state.atoms[a].fixed_directions() for a in down - sub}
                after = state.cut_as_free(sub)
                for a in down - sub:
                    if (
                        Direction.BOTTOM in after.atoms[a].fixed_directions()
                        and Direction.BOTTOM not in before[a]
                    ):
                        bottom_source[a] = idx
                state = after
            if found_no_slack:
                break
        assert found_no_slack, "未能复现该 charging 方案的失败案例；FINDINGS 中的负结果需复核"


class TestInvarianceUnderTieBreaks:
    """[数值证据·强] g 与 n33 在所有 tie-break 分支间恒定（FINDINGS F12）。

    若成立，则 Σb_i 也是 molecule 的不变量，命题 P 遂成为两个不变量之间的
    不等式，而非关于执行过程的摊还界——这改变了证明的目标形状。

    同时这是拟阵结构的**必要条件**（所有极大可行解等基数）。
    """

    @pytest.mark.parametrize("n", list(range(6, 14)))
    def test_g_and_n33_constant_across_branches(self, n):
        from molecule_cut.exhaustive import enumerate_tiebreak_records

        rng = random.Random(31337 + n)
        multi = 0
        for _ in range(25):
            mol = random_toy1_full(n, rng.randint(max(1, n - 2), n + 3), rng)
            if mol is None:
                continue
            try:
                recs = enumerate_tiebreak_records(mol, cap=3000)
            except RuntimeError:
                continue
            gs, n33s = set(), set()
            for _, rec in recs:
                if rec.failed:
                    continue
                c = rec.op_counts()
                gs.add(c.get(OpKind.CUT33, 0) + c.get(OpKind.CUT343, 0))
                n33s.add(rec.n33)
            if not gs:
                continue
            if len([r for _, r in recs if r.failed is None]) > 1:
                multi += 1
            assert len(gs) == 1, f"n={n}: g 在分支间取值 {sorted(gs)}"
            assert len(n33s) == 1, f"n={n}: n33 在分支间取值 {sorted(n33s)}"
        # 守卫（LESSONS 13）：确认确实比较了多分支实例
        assert multi > 0, f"n={n}: 没有多分支实例，不变性测试形同空转"
