"""穷举 tie-breaking，计算 a_n；以及全搜合法切割序列，计算 v_d(M)。

SPEC §1.6 的三个量：
  a_n     = min over M, min over tie-breaks, #{33}（Def 1.4 算法的最坏值）
  v_d(M)  = max over legal cuts, (#{33} - 10d·#{4})（molecule 的内在量）

tie-breaking 位点（SPEC §3.6）：step1（多个最低 MU atom）、step2（多条跨层
bond，Toy I 中无）、step3（多对相邻 deg 3）、step4（多个 {343}）、step6
（多个 deg 3）。
"""

from __future__ import annotations

from dataclasses import dataclass

from .algorithm import CutRecord, TieBreaker, run_toy1_algorithm
from .elementary import ElemType, classify, is_elementary
from .molecule import Molecule


class ScriptedTieBreaker(TieBreaker):
    """按预先给定的索引序列做选择；索引用尽时回退到字典序最小。

    与 `enumerate_tiebreak_records` 配合，用迭代加深的方式穷举所有分支。
    """

    def __init__(self, script: tuple[int, ...]) -> None:
        self.script = script
        self.pos = 0
        self.trace: list[int] = []  # 实际遇到的每个选择点的候选数

    def pick(self, site: str, options):
        opts = sorted(options, key=lambda s: sorted(s) if isinstance(s, set) else s)
        self.trace.append(len(opts))
        if len(opts) == 1:
            return opts[0]
        idx = self.script[self.pos] if self.pos < len(self.script) else 0
        self.pos += 1
        return opts[min(idx, len(opts) - 1)]


def enumerate_tiebreak_records(mol: Molecule, cap: int = 200_000):
    """穷举所有 tie-breaking 分支，产出 (script, CutRecord)。

    用 DFS：先跑一次拿到第一个选择点的候选数，再对每个分支递归。
    cap 限制探索的分支总数，超出时抛错（避免静默截断——SPEC 纪律）。
    """
    results: list[tuple[tuple[int, ...], CutRecord]] = []
    stack: list[tuple[int, ...]] = [()]
    seen: set[tuple[int, ...]] = set()
    explored = 0

    while stack:
        script = stack.pop()
        if script in seen:
            continue
        seen.add(script)
        explored += 1
        if explored > cap:
            raise RuntimeError(f"tie-break exploration exceeded cap={cap}")

        tb = ScriptedTieBreaker(script)
        rec = run_toy1_algorithm(mol, choose=tb)
        results.append((script, rec))

        # 在 script 用尽之后的第一个多选点上分叉
        multi_positions = [i for i, k in enumerate(tb.trace) if k > 1]
        for pos_idx, pos in enumerate(multi_positions):
            if pos_idx < len(script):
                continue
            n_opts = tb.trace[pos]
            for alt in range(1, n_opts):
                stack.append(script + (alt,))
            break

    return results


@dataclass
class AnResult:
    n33_min: int
    n33_max: int
    n4_values: set[int]
    n_branches: int
    worst_script: tuple[int, ...]


def a_of_molecule(mol: Molecule) -> AnResult:
    """对单个 molecule 穷举 tie-breaking，返回 #{33} 的最小/最大值。"""
    recs = enumerate_tiebreak_records(mol)
    # a_n 是 min over tie-breaks；静默丢弃失败分支只会把 min 抬高。Prop 1.6 保证
    # 合法 Toy I 输入下每个分支都完成，故任一 failed 都是实现 bug 或非 Toy I 输入，
    # 必须立刻暴露而不是被过滤掉。
    bad = [(s, r) for s, r in recs if r.failed is not None]
    if bad:
        script, rec = bad[0]
        raise RuntimeError(
            f"tie-break branch {script} failed for {mol}: {rec.failed} "
            f"({len(bad)} of {len(recs)} branches failed)"
        )
    ok = recs
    n33s = [r.n33 for _, r in ok]
    lo = min(n33s)
    worst = next(s for s, r in ok if r.n33 == lo)
    return AnResult(
        n33_min=lo,
        n33_max=max(n33s),
        n4_values={r.n4 for _, r in ok},
        n_branches=len(ok),
        worst_script=worst,
    )


# ---------- v_d(M)：全搜合法切割序列 ----------


def _legal_moves(mol: Molecule) -> list[frozenset[int]]:
    """当前 molecule 上所有能切出 elementary 组件的子集。

    候选限制为 1 或 2 个 atom（elementary molecule 至多 2 个 atom），
    以及 {343} 三元组（可再切成 {3}+{33}）。
    """
    moves: list[frozenset[int]] = []
    ids = sorted(mol.atoms)

    for a in ids:
        if is_elementary(mol, {a}):
            moves.append(frozenset({a}))

    for b in mol.bonds:
        pair = {b.parent, b.child}
        if is_elementary(mol, pair):
            moves.append(frozenset(pair))

    # {343}：作为一步产出 {3} + {33}
    for a in ids:
        if mol.degree(a) != 4:
            continue
        deg3 = sorted(nb for nb in mol.neighbors(a) if mol.degree(nb) == 3)
        for i in range(len(deg3)):
            for j in range(i + 1, len(deg3)):
                triple = {deg3[i], a, deg3[j]}
                if classify(mol, triple) is ElemType.E343:
                    moves.append(frozenset(triple))

    return list(dict.fromkeys(moves))


def v_d(mol: Molecule, d: int = 3, memo: dict | None = None) -> int | None:
    """molecule 的内在最优切割值（SPEC §1.6）。

    默认走 `fast_vd.fast_v_d`（位掩码 DP）。`memo` 参数保留以兼容旧调用，
    但在默认路径上被忽略。需要可逐步调试的对象级实现时用 `v_d_reference`。
    """
    if memo is not None:
        return v_d_reference(mol, d=d, memo=memo)
    from .fast_vd import fast_v_d

    return fast_v_d(mol, d=d)


def v_d_reference(mol: Molecule, d: int = 3, memo: dict | None = None) -> int | None:
    """慢速参考实现：每步构造 `Molecule` 子图。仅用于差分测试。

    返回 max over legal cutting sequences of (#{33} - 10d·#{4})；
    若不存在把 M 完全分解为 elementary 组件的合法序列，返回 None。
    """
    memo = {} if memo is None else memo
    key = _state_key(mol)
    if key in memo:
        return memo[key]

    if not mol.atoms:
        memo[key] = 0
        return 0

    best: int | None = None
    for mv in _legal_moves(mol):
        t = classify(mol, set(mv))
        if t is ElemType.E33:
            gain = 1
        elif t is ElemType.E4:
            gain = -10 * d
        elif t is ElemType.E343:
            gain = 1  # {3} + {33}
        else:
            gain = 0
        sub = v_d_reference(mol.cut_as_free(set(mv)), d, memo)
        if sub is None:
            continue
        val = gain + sub
        if best is None or val > best:
            best = val

    memo[key] = best
    return best


def w(mol: Molecule, memo: dict | None = None) -> int | None:
    """内在好分量数 `W(M) = max_C #{33}(C)`。

    `W` 与 `v_d + 10d` 只有在某个 W 最优序列恰有一个 `{4}` 时才相同；
    此函数故意单独求解，避免把该条件当作默认事实。
    """
    if memo is not None:
        return w_reference(mol, memo=memo)
    from .fast_vd import fast_w

    return fast_w(mol)


def w_reference(mol: Molecule, memo: dict | None = None) -> int | None:
    """慢速参考实现：直接枚举合法切割序列来计算 `W`。"""
    memo = {} if memo is None else memo
    key = _state_key(mol)
    if key in memo:
        return memo[key]
    if not mol.atoms:
        memo[key] = 0
        return 0

    best: int | None = None
    for mv in _legal_moves(mol):
        t = classify(mol, set(mv))
        gain = 1 if t in {ElemType.E33, ElemType.E343} else 0
        sub = w_reference(mol.cut_as_free(set(mv)), memo)
        if sub is None:
            continue
        val = gain + sub
        if best is None or val > best:
            best = val
    memo[key] = best
    return best


def _state_key(mol: Molecule) -> tuple:
    """molecule 状态的哈希键（用于 v_d 的记忆化）。

    注意：这里用的是带标号的键，不是同构不变量——同构状态会被重复计算。
    正确但不最优；n 较小时可接受。
    """
    atoms = tuple(
        (
            aid,
            mol.atoms[aid].layer.value,
            tuple(sorted((k.value, dd.value) for k, dd in mol.atoms[aid].ends)),
        )
        for aid in sorted(mol.atoms)
    )
    bonds = tuple(sorted((b.parent, b.child) for b in mol.bonds))
    return (atoms, bonds)
