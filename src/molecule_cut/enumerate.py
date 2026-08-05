"""Toy model I 的受限枚举与小规模全域覆盖枚举（SPEC §1.5）。

Toy I 的结构约束（Def 1.2 + Def 1.1）：
  * MU、MD 各为树，无 fixed end
  * 每个 MU atom 恰有一条 bond 连到 MD（#2conn = 0 且 X = MU）
  * MD 中无 atom 是 MU 中 atom 的 parent（跨层 bond 方向恒为 U→D）
  * 每个 atom 至多 2 parent、至多 2 child

受限枚举器的槽位推论：
  * MU atom：1 个 child 槽被跨层 bond 占用 → MU 内至多 1 child
    但它仍可有 2 个 parent，所以完整定义域中的 MU 可以是倒置二叉树，
    **不必是链**。`enumerate_toy1*` 为快速实验固定 MU 为链；
    `enumerate_toy1_full_labeled` 才枚举两层的全部合法定向。
  * MD atom：至多 2 parent，其中跨层 bond 占用 parent 槽
    ⟹ MD 内 parent 数 + 跨层入度 ≤ 2。

因此受限 Toy I 的枚举 = （MD 的单父树形状）×（|MU| 条跨层 bond 到 MD 的分配）。
同构去重：MU 是链故无对称性，只需对 MD 的树同构 + 跨层分配一起做
canonical form。
"""

from __future__ import annotations

from collections.abc import Iterator
from itertools import product

from .builders import build
from .molecule import Molecule


def rooted_binary_trees(n: int) -> Iterator[tuple[int, ...]]:
    """生成 n 个节点的有根树，每节点至多 2 child，节点 0 为根。

    表示为 parent 数组 `par`，par[i] = i 的 parent（par[0] = -1），
    且强制 par[i] < i（保证无环、根为 0）。

    为避免同构重复，对同一 parent 的多个 child 施加"标签递增"约束不足以
    完全去重，故这里生成所有 labeled 形态，由 canonical_form 统一去重。
    """
    if n <= 0:
        return
    if n == 1:
        yield (-1,)
        return

    def rec(i: int, par: list[int], child_count: list[int]) -> Iterator[tuple[int, ...]]:
        if i == n:
            yield tuple(par)
            return
        for p in range(i):
            if child_count[p] >= 2:
                continue
            par.append(p)
            child_count[p] += 1
            child_count.append(0)
            yield from rec(i + 1, par, child_count)
            child_count.pop()
            child_count[p] -= 1
            par.pop()

    yield from rec(1, [-1], [0])


def labelled_trees(n: int) -> Iterator[tuple[tuple[int, int], ...]]:
    """生成 n 个标号节点的全部无向树（Prüfer 序列，每棵一次）。"""
    if n <= 0:
        return
    if n == 1:
        yield ()
        return
    if n == 2:
        yield ((0, 1),)
        return

    for sequence in product(range(n), repeat=n - 2):
        degree = [1] * n
        for vertex in sequence:
            degree[vertex] += 1
        edges: list[tuple[int, int]] = []
        for vertex in sequence:
            leaf = next(index for index, value in enumerate(degree) if value == 1)
            edges.append((leaf, vertex))
            degree[leaf] -= 1
            degree[vertex] -= 1
        rest = [index for index, value in enumerate(degree) if value == 1]
        edges.append((rest[0], rest[1]))
        yield tuple(edges)


def oriented_trees(
    n: int, *, max_parents: int = 2, max_children: int = 2
) -> Iterator[tuple[tuple[tuple[int, int], ...], tuple[int, ...], tuple[int, ...]]]:
    """全部标号树及其满足 parent/child 槽位上限的定向。

    不假设单父根树。对 MU，调用方设 `max_children=1`，为每个 atom
    留出一条跨层 child bond；这仍允许倒置二叉树中的双父 atom。
    """
    for edges in labelled_trees(n):
        for bits in product((0, 1), repeat=len(edges)):
            oriented = tuple(
                (left, right) if bit == 0 else (right, left)
                for (left, right), bit in zip(edges, bits, strict=True)
            )
            n_parents = [0] * n
            n_children = [0] * n
            for parent, child in oriented:
                n_children[parent] += 1
                n_parents[child] += 1
            if max(n_parents, default=0) > max_parents or max(n_children, default=0) > max_children:
                continue
            yield oriented, tuple(n_parents), tuple(n_children)


def canonical_form(mol: Molecule, n_up: int, n_down: int) -> tuple:
    """Toy I molecule 的规范形式，用于同构去重。

    MU 是链（0 为顶、n_up-1 为最低），无对称性，故 MU 的编号固定。
    MD 的编号可任意置换，因此用"MD 树的 AHU 规范编码 + 每个 MD 节点收到
    的跨层 bond 来自哪个 MU 位置"联合构造不变量。

    实现：对 MD 树做自底向上的规范编码，节点的标签包含其跨层入边集合。
    """
    up_ids = list(range(n_up))
    down_ids = [n_up + j for j in range(n_down)]

    # MD 内部的 parent/child 关系
    md_children: dict[int, list[int]] = {d: [] for d in down_ids}
    md_parent: dict[int, int | None] = {d: None for d in down_ids}
    cross_from: dict[int, list[int]] = {d: [] for d in down_ids}

    for b in mol.bonds:
        if b.parent in up_ids and b.child in down_ids:
            cross_from[b.child].append(up_ids.index(b.parent))
        elif b.parent in down_ids and b.child in down_ids:
            md_children[b.parent].append(b.child)
            md_parent[b.child] = b.parent

    roots = [d for d in down_ids if md_parent[d] is None]

    def encode(node: int) -> tuple:
        """AHU 风格递归编码：(跨层入边签名, 排序后的子树编码)。"""
        sig = tuple(sorted(cross_from[node]))
        subs = sorted(encode(c) for c in md_children[node])
        return (sig, tuple(subs))

    return (n_up, tuple(sorted(encode(r) for r in roots)))


def enumerate_toy1(n: int) -> Iterator[Molecule]:
    """枚举 |MU| = n 的受限非同构 Toy I 子类。

    MU 为 n 节点链；MD 为 n 节点树（Def 1.1(3) 未强制 |MD| = |MU|，
    但 Toy I 中每个 MU atom 恰有一条跨层 bond，若 |MD| < n 则某些 MD atom
    需承接多条跨层 bond。这里先取 |MD| = n 的主情形，|MD| != n 由
    enumerate_toy1_general 处理）。
    """
    yield from enumerate_toy1_general(n, n)


def enumerate_toy1_general(n_up: int, n_down: int) -> Iterator[Molecule]:
    """受限子类：MU 为链、MD 为单父有根树。

    这是快速实验用的 canonical 子类，**不是** Toy I 完整定义域；完整标号
    枚举见 `enumerate_toy1_full_labeled`。
    """
    up_edges = [(i, i + 1) for i in range(n_up - 1)]
    seen: set[tuple] = set()

    for par in rooted_binary_trees(n_down):
        down_edges = [(par[i], i) for i in range(1, n_down)]
        # 每个 MD 节点在 MD 内已用掉的 parent 槽
        md_indeg = [0] * n_down
        for i in range(1, n_down):
            md_indeg[i] = 1
        # 剩余 parent 槽（至多 2）
        capacity = [2 - md_indeg[j] for j in range(n_down)]

        # 把 n_up 条跨层 bond 分配到 MD 节点，受 capacity 限制
        for assign in product(range(n_down), repeat=n_up):
            used = [0] * n_down
            ok = True
            for tgt in assign:
                used[tgt] += 1
                if used[tgt] > capacity[tgt]:
                    ok = False
                    break
            if not ok:
                continue

            cross = [(u, assign[u]) for u in range(n_up)]
            try:
                mol = build(n_up, n_down, up_edges, down_edges, cross)
            except ValueError:
                continue  # 槽位冲突

            cf = canonical_form(mol, n_up, n_down)
            if cf in seen:
                continue
            seen.add(cf)
            yield mol


def enumerate_toy1_full_labeled(n_up: int, n_down: int) -> Iterator[Molecule]:
    """完整 Toy I 的标号枚举，不做同构去重。

    MU/MD 都遍历全部无向树和合法定向。MU 的层内 child 数限制为 1，
    因为每个 MU atom 还要有一条跨层 child bond；它仍可有两个 layer-internal
    parents，故不被错误地限制为链。此生成器只适用于很小规模的覆盖审计。
    """
    for up_edges, _, _ in oriented_trees(n_up, max_children=1):
        for down_edges, down_parents, _ in oriented_trees(n_down):
            capacity = [2 - count for count in down_parents]
            for assignment in product(range(n_down), repeat=n_up):
                used = [0] * n_down
                for target in assignment:
                    used[target] += 1
                if any(used[index] > capacity[index] for index in range(n_down)):
                    continue
                yield build(
                    n_up,
                    n_down,
                    list(up_edges),
                    list(down_edges),
                    list(enumerate(assignment)),
                )


def count_toy1(n: int) -> int:
    return sum(1 for _ in enumerate_toy1(n))
