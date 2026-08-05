"""显式构造的 molecule 族。

核心族是 `transpose_family(k)`（n = 3k+1），命名来自其跨层 bond 分配是
3×k 网格的转置映射。另有两个变体覆盖其余剩余类，合起来给出 FINDINGS F9 的
上界 a_n ≤ ⌈(n−1)/3⌉ 对所有 n ≥ 2：

    n ≡ 1 (mod 3)   transpose_family(k)          n = 3k+1
    n ≡ 2 (mod 3)   padded_family(k)             n = 3k+2
    n ≡ 0 (mod 3)   trimmed_family(k)            n = 3k
"""

from __future__ import annotations

from .builders import build
from .molecule import Molecule


def transpose_cross_target(i: int, k: int) -> int:
    """MU atom u_i 的跨层 bond 落在哪个 MD atom（局部编号）。

    对 0 ≤ i < 3k，写 i = q·k + r（q ∈ {0,1,2}, 0 ≤ r < k），映到 3r + q；
    这是 3×k 网格的转置，是 [0,3k) 上的双射。额外的 u_{3k} 落在 MD_0。

    效果：MD 路径上的第 r 个三元块 {3r, 3r+1, 3r+2} 恰好收到来自
    u_r、u_{k+r}、u_{2k+r} 的三条跨层 bond——即 MU 链三个「三分之一」各一条。
    """
    if not 0 <= i <= 3 * k:
        raise ValueError(f"i={i} out of range [0, {3 * k}]")
    if i == 3 * k:
        return 0
    q, r = divmod(i, k)
    return 3 * r + q


def transpose_family(k: int) -> Molecule:
    """构造 |MU| = 3k+1、|MD| = 3k 的 toy model I 实例。

    结构：
      * MU：3k+1 个 atom 的链 u_0 → u_1 → … → u_{3k}（u_{3k} 最低，最先被切）
      * MD：3k 个 atom 的路径 d_0 → d_1 → … → d_{3k-1}
      * 跨层：u_i → d_{transpose_cross_target(i, k)}

    在 Def 1.4 算法下的执行（已对 k = 1..12 计算验证，见 tests）：
      阶段 1  k+1 次操作 (a)，依次切 u_{3k}, u_{3k-1}, …, u_{2k}
      阶段 2  k 次 {343}，依次切 MD 三元块 [0,1,2], [3,4,5], …（级联触发）
      阶段 3  2k 次 CUT_U_ONLY，切剩余的 u_{2k-1}, …, u_0
    产出 #{33} = k，操作剖面恒为 (A,B,C,D,cutU) = (k+1, 0, 0, k, 2k)。
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    n = 3 * k + 1
    md = 3 * k
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = [(j, j + 1) for j in range(md - 1)]
    cross = [(i, transpose_cross_target(i, k)) for i in range(n)]
    return build(n, md, up_edges, down_edges, cross)


def padded_family(k: int) -> Molecule:
    """n ≡ 2 (mod 3) 的变体：n = 3k+2，|MD| = 3k+1。

    在 `transpose_family(k)` 的基础上，于 MU 链**顶端**追加 1 个 atom、
    MD 路径**尾端**追加 1 个 atom，并把新 MU atom 接到新 MD atom 上。
    顶端的 MU atom 最后才被切，故不干扰阶段 1/2 的级联。

    产出 #{33} = k+1 = ⌈(n−1)/3⌉（已对 k = 1..8 计算验证）。
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    n = 3 * k + 2
    md = 3 * k + 1
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = [(j, j + 1) for j in range(md - 1)]
    cross = [(0, 3 * k)] + [(i, transpose_cross_target(i - 1, k)) for i in range(1, n)]
    return build(n, md, up_edges, down_edges, cross)


def trimmed_family(k: int) -> Molecule:
    """n ≡ 0 (mod 3) 的变体：n = 3k，|MD| = 3k−1。

    去掉 `transpose_family(k)` 中额外的 u_{3k}，并把 MD 路径缩短一个；
    原本落在被删 MD atom 上的跨层 bond 改接到 d_0（使 d_0 仍收 2 条）。

    产出 #{33} = k = ⌈(n−1)/3⌉（已对 k = 2..8 计算验证）。
    """
    if k < 2:
        raise ValueError("k must be >= 2")
    n = 3 * k
    md = 3 * k - 1
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = [(j, j + 1) for j in range(md - 1)]
    cross = []
    for i in range(n):
        t = transpose_cross_target(i, k)
        cross.append((i, t if t < md else 0))
    return build(n, md, up_edges, down_edges, cross)


def family_for(n: int) -> Molecule:
    """按 n 的剩余类返回相应族的实例。n >= 4。"""
    if n < 4:
        raise ValueError("n must be >= 4")
    q, s = divmod(n, 3)
    if s == 1:
        return transpose_family(q)
    if s == 2:
        return padded_family(q)
    return trimmed_family(q)
