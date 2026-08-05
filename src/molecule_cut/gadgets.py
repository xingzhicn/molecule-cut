"""B2/X3C spike 用的手造 molecule 族。

这些不是定理构造，只用于数值探测 hardness gadget 是否可能工作。
全部经 `build` 进入合法简化两层（MU/MD 树 + 槽位约束）。
"""

from __future__ import annotations

from .builders import build
from .molecule import Molecule


def hubs_2conn(m: int) -> Molecule:
    """m 个 #2conn 叶挂在 hub 路径上；每个叶双跨到 MD 路径上的私有边。

    |U|=|D|=2m。用于观察 W 与 rho / #2conn 的关系。
    """
    if m < 1:
        raise ValueError("m >= 1")
    u, d = 2 * m, 2 * m
    up_edges = [(i, i + 1) for i in range(m - 1)] + [(i, m + i) for i in range(m)]
    down_edges = [(i, i + 1) for i in range(d - 1)]
    cross: list[tuple[int, int]] = []
    for i in range(m):
        cross.append((m + i, 2 * i))
        cross.append((m + i, 2 * i + 1))
    return build(u, d, up_edges, down_edges, cross)


def element_spine_md(n_el: int) -> tuple[int, list[tuple[int, int]]]:
    """MD：n_el 个元素叶 + n_el 个 spine 节点；每个元素恰 1 个 MD 父 → 跨层 cap≤1。

    返回 (n_down, down_edges)。元素 id = 0..n_el-1，spine id = n_el..2*n_el-1。
    """
    if n_el < 1:
        raise ValueError("n_el >= 1")
    s0 = n_el
    down_edges = [(s0 + i, s0 + i + 1) for i in range(n_el - 1)]
    down_edges += [(s0 + i, i) for i in range(n_el)]
    return 2 * n_el, down_edges


def three_set_gadget() -> Molecule:
    """单个 3-集 gadget：hub + 2conn 叶 + 1conn 叶，盖住 3 个 cap-1 元素。"""
    n_down, de = element_spine_md(3)
    # U: hub=0, L=1 (2conn), R=2 (1conn)
    return build(3, n_down, [(0, 1), (0, 2)], de, [(1, 0), (1, 1), (2, 2)])


def two_set_geometry(
    cross: list[tuple[int, int]],
) -> Molecule:
    """固定两-集 MU/MD 骨架，只换跨层（用于 YES / half / scramble 对照）。

    U: root=0, hA=1, hB=2, LA=3, RA=4, LB=5, RB=6。
    元素 0..5 为 MD 叶（cap 1）。
    """
    n_down, de = element_spine_md(6)
    ue = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]
    return build(7, n_down, ue, de, cross)


def two_sets_yes() -> Molecule:
    """不相交两 3-集：{0,1,2} 与 {3,4,5}。"""
    return two_set_geometry(
        [(3, 0), (3, 1), (4, 2), (5, 3), (5, 4), (6, 5)]
    )


def two_sets_half() -> Molecule:
    """同一骨架只激活第一集。"""
    return two_set_geometry([(3, 0), (3, 1), (4, 2)])


def two_sets_scramble() -> Molecule:
    """同数量跨层、打乱元素归属（非两个干净 3-集）。"""
    return two_set_geometry(
        [(3, 0), (3, 3), (4, 1), (5, 2), (5, 4), (6, 5)]
    )


def two_sets_overlap_cross() -> list[tuple[int, int]]:
    """故意 overlap 的跨层表（期望 build 因元素 cap 失败）。"""
    return [(3, 0), (3, 1), (4, 2), (5, 1), (5, 2), (6, 3)]
