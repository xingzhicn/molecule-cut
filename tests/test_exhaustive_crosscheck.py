"""tie-breaking 穷举的对拍测试。

`enumerate_tiebreak_records` 用迭代加深 DFS 探索选择树，逻辑不显然。
这里用一个"明显正确但低效"的递归枚举器做差分测试：它在每个选择点直接
分叉，不依赖 script 编码。两者产出的 #{33} 多重集必须完全一致。
"""

import pytest

from molecule_cut.algorithm import TieBreaker, run_toy1_algorithm
from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.exhaustive import enumerate_tiebreak_records


class ForcedChoice(TieBreaker):
    """在第 k 个多选点强制选第 idx 个候选，其余点选 0。用于朴素枚举。"""

    def __init__(self, prefix: tuple[int, ...]) -> None:
        self.prefix = prefix
        self.seen_multi = 0
        self.branch_points: list[int] = []

    def pick(self, site, options):
        opts = sorted(options, key=lambda s: sorted(s) if isinstance(s, set) else s)
        if len(opts) == 1:
            return opts[0]
        k = self.seen_multi
        self.seen_multi += 1
        self.branch_points.append(len(opts))
        if k < len(self.prefix):
            return opts[self.prefix[k]]
        return opts[0]


def naive_all_n33(mol) -> list[int]:
    """朴素递归：逐层展开所有选择组合，返回所有分支的 #{33} 列表。"""
    out: list[int] = []
    stack: list[tuple[int, ...]] = [()]
    done: set[tuple[int, ...]] = set()

    while stack:
        prefix = stack.pop()
        if prefix in done:
            continue
        done.add(prefix)

        tb = ForcedChoice(prefix)
        rec = run_toy1_algorithm(mol, choose=tb)
        if rec.failed is None:
            out.append(rec.n33)

        # 若 prefix 已覆盖所有多选点，则这是一个叶子；否则在下一个点分叉
        if len(prefix) < len(tb.branch_points):
            k = len(prefix)
            for alt in range(tb.branch_points[k]):
                stack.append(prefix + (alt,))

    return out


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_tiebreak_enumeration_matches_naive(n):
    """两套枚举器给出相同的 #{33} 取值集合与最小值。

    **覆盖漏洞记录（2026-07-29 审查发现）**：本测试原先只参数化到 n=4，而
    n ≤ 4 的每个 Toy I 实例都**只有一条** tie-break 分支——两个枚举器各自
    只产出一个结果就比对通过，**分支逻辑从未被真正验证**。多分支实例最早在
    n=5 出现，故必须覆盖到 n=5；`test_branching_is_actually_exercised`
    进一步断言该覆盖不会再退化成空转。
    """
    for idx, mol in enumerate(enumerate_toy1(n)):
        fast = [r.n33 for _, r in enumerate_tiebreak_records(mol) if r.failed is None]
        slow = naive_all_n33(mol)

        assert fast, f"n={n} mol#{idx}: fast enumerator produced nothing"
        assert slow, f"n={n} mol#{idx}: naive enumerator produced nothing"
        assert min(fast) == min(slow), (
            f"n={n} mol#{idx}: min mismatch fast={min(fast)} slow={min(slow)}"
        )
        assert set(fast) == set(slow), (
            f"n={n} mol#{idx}: value set mismatch fast={sorted(set(fast))} slow={sorted(set(slow))}"
        )


def test_branching_is_actually_exercised():
    """守卫：确保上面的对拍确实比较了多分支实例，而不是空转通过。

    若某次重构使多分支实例消失（例如 tie-break 位点被误删），本测试会失败，
    从而暴露「对拍看似通过但什么都没验证」的情况。
    """
    multi = 0
    max_branches = 0
    for mol in enumerate_toy1(5):
        k = len([r for _, r in enumerate_tiebreak_records(mol) if r.failed is None])
        max_branches = max(max_branches, k)
        if k > 1:
            multi += 1
    assert multi > 0, "n=5 上没有任何多分支实例，对拍测试形同虚设"
    assert max_branches >= 2, f"最大分支数只有 {max_branches}"


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_n4_is_one_across_all_tiebreaks(n):
    """Prop 1.6 的 #{4} = 1 在所有 tie-breaking 分支下都成立。"""
    for mol in enumerate_toy1(n):
        for _, rec in enumerate_tiebreak_records(mol):
            if rec.failed is None:
                assert rec.n4 == 1


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_prop16_bound_across_all_tiebreaks(n):
    """M1 正确性 Gate：所有 molecule、所有 tie-break 下 #{33} ≥ (n-1)/5。"""
    for mol in enumerate_toy1(n):
        for _, rec in enumerate_tiebreak_records(mol):
            if rec.failed is None:
                assert rec.n33 >= (n - 1) / 5
