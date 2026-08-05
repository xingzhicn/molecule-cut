"""Def 1.4：toy model I 的切割算法，以及切割过程的记录。

严格实现 SPEC.md §3.6 的六个步骤，并显式暴露其中四处 tie-breaking 自由度
（这正是 a_n 定义里 "min over tie-breaks" 的含义）。

同时记录 SPEC §3.7 的 λ 势函数与每步的操作类型，用于复现 Prop 1.6 的计数
论证（这是最强的正确性校验）。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from .elementary import (
    ElemType,
    classify,
    find_343_triples,
    find_adjacent_deg3_pairs,
    is_proper,
    lambda_potential,
    violates_monotonicity,
)
from .molecule import Layer, Molecule


class OpKind(Enum):
    """SPEC §3.7 表中的操作类型。"""

    A = "(a)"  # 切 n as free，m 存在且为 deg 4 —— λ += 1，最多断 1 条跨层 bond
    B = "(b)"  # 切 {n,m} as free，m 为 deg 3 —— λ 不变，最多断 2 条
    CUT33 = "cut{33}"  # λ -= 1，最多断 3 条
    CUT343 = "cut{343}"  # λ -= 1，最多断 4 条
    CUT3 = "cut{3}"  # λ 不变
    CUT_U_ONLY = "cut_U"  # 切一个 MD 邻居已被移除的 MU atom —— λ 不变，断 0 条
    # 注：Def 1.4 步骤 (2) 的措辞 "otherwise just cut n as free" 同时涵盖
    # (a) 与 CUT_U_ONLY 两种情形，但二者对 λ 的作用不同：(a) 把某个 deg 4
    # 的 m 变成 deg 3（λ += 1），而 CUT_U_ONLY 完全不触及 MD（λ 不变）。
    # 讲义 Prop 1.6 证明 part 1 只列举了前两种情形，未单列 CUT_U_ONLY；
    # 恒等式 (1.4) 要成立，必须把 (a) 严格限定为「m 存在且为 deg 4」。
    # 实测：Toy I 中 CUT_U_ONLY 确实出现（当 n 的 MD 邻居已被此前的
    # {33}/{343} 带走时），产出 {2} molecule（normal，计分 0）。


@dataclass
class Step:
    """一步切割的完整记录。"""

    op: OpKind
    subset: frozenset[int]
    elem_types: tuple[ElemType, ...]  # 该步产出的 elementary 组件类型
    lam_before: int
    lam_after: int
    cross_bonds_broken: int


@dataclass
class CutRecord:
    """一次完整切割的结果。"""

    steps: list[Step] = field(default_factory=list)
    n33: int = 0
    n4: int = 0
    n3: int = 0
    n2: int = 0
    failed: str | None = None  # 非 None 表示中途违反不变量

    def score(self, d: int) -> int:
        """SPEC §1.2：#{33} - 10d * #{4}。"""
        return self.n33 - 10 * d * self.n4

    def op_counts(self) -> dict[OpKind, int]:
        out: dict[OpKind, int] = {}
        for s in self.steps:
            out[s.op] = out.get(s.op, 0) + 1
        return out


def _count_cross_broken(mol: Molecule, subset: set[int]) -> int:
    """subset 被切走时打断的跨层 bond 数（用于校验 SPEC §3.7 的上界表）。"""
    n = 0
    for b in mol.cross_bonds():
        p_in, c_in = b.parent in subset, b.child in subset
        if p_in != c_in or (p_in and c_in):
            n += 1
    return n


def _record(
    rec: CutRecord, mol: Molecule, subset: set[int], op: OpKind, pieces: list[set[int]]
) -> Molecule:
    """执行一步切割：分类各产出组件、更新计数、返回切割后的 molecule。"""
    lam_before = lambda_potential(mol)
    types = tuple(classify(mol, p) for p in pieces)
    broken = _count_cross_broken(mol, subset)

    new_mol = mol.cut_as_free(subset)
    lam_after = lambda_potential(new_mol)

    for t in types:
        if t is ElemType.E33:
            rec.n33 += 1
        elif t is ElemType.E4:
            rec.n4 += 1
        elif t is ElemType.E3:
            rec.n3 += 1
        elif t is ElemType.E2:
            rec.n2 += 1
        elif t is ElemType.E343:
            # Def 1.4 (4)：{343} 再切成 {3} + {33}，计一个 {33}（见 (1.3)）
            rec.n33 += 1
            rec.n3 += 1
        else:
            rec.failed = f"non-elementary component {t} at subset {sorted(subset)}"

    rec.steps.append(Step(op, frozenset(subset), types, lam_before, lam_after, broken))
    return new_mol


def _lowest_uncut_up(mol: Molecule) -> list[int]:
    """Def 1.4 步骤 (1)：MU 中「最低」的未切 atom —— 在 MU 内无 child。

    返回所有候选（tie-breaking 自由度之一；Toy I 初始时应恰有一个）。
    """
    up = set(mol.layer_atoms(Layer.UP))
    return [a for a in sorted(up) if not (set(mol.children(a)) & up)]


def _cleanup_md(mol: Molecule, rec: CutRecord, choose: TieBreaker) -> Molecule:
    """Def 1.4 步骤 (3)-(4)：反复取走 {33} 与 {343}，直到 MD 变 proper。"""
    while True:
        pairs = find_adjacent_deg3_pairs(mol)
        if pairs:
            sel = choose.pick("step3", pairs)
            mol = _record(rec, mol, sel, OpKind.CUT33, [sel])
            continue
        triples = find_343_triples(mol)
        if triples:
            sel = choose.pick("step4", triples)
            mol = _record(rec, mol, sel, OpKind.CUT343, [sel])
            continue
        return mol


class TieBreaker:
    """封装 SPEC §3.6 列出的四处 tie-breaking 选择。

    默认实现取字典序最小（确定性）。枚举器会用穷举版本覆盖所有选择。
    """

    def pick(self, site: str, options: list[set[int]] | list[int]):
        if not options:
            raise ValueError(f"no options at {site}")
        if isinstance(options[0], set):
            return min(options, key=lambda s: sorted(s))  # type: ignore[arg-type]
        return min(options)  # type: ignore[type-var]


def run_toy1_algorithm(
    mol: Molecule, choose: TieBreaker | None = None, check: bool = True
) -> CutRecord:
    """对 toy model I 运行 Def 1.4 的算法。

    Args:
        mol: 输入 molecule（应为 toy model I）
        choose: tie-breaking 策略
        check: 是否在每步后断言 monotonicity（SPEC §3.4）
    """
    choose = choose or TieBreaker()
    rec = CutRecord()
    cur = mol.copy()

    guard = 0
    max_steps = 100 * (len(mol.atoms) + 1)

    while cur.layer_atoms(Layer.UP):
        guard += 1
        if guard > max_steps:
            rec.failed = "step limit exceeded (possible infinite loop)"
            return rec

        # 步骤 (1)/(5)：选一个最低的未切 MU atom
        candidates = _lowest_uncut_up(cur)
        if not candidates:
            rec.failed = "no lowest atom in MU but MU non-empty"
            return rec
        n = choose.pick("step1", candidates)

        # 步骤 (2)：看 n 在 MD 的邻居
        down = set(cur.layer_atoms(Layer.DOWN))
        md_nbrs = [nb for nb in cur.children(n) if nb in down]

        if not md_nbrs:
            cur = _record(rec, cur, {n}, OpKind.CUT_U_ONLY, [{n}])
        else:
            m = choose.pick("step2", md_nbrs)  # Toy I 中唯一，无自由度
            if cur.degree(m) == 3:
                # 操作 (b)：切 {n, m} as free
                cur = _record(rec, cur, {n, m}, OpKind.B, [{n, m}])
            else:
                # 操作 (a)：只切 n as free（m 为 deg 4，将变 deg 3）
                cur = _record(rec, cur, {n}, OpKind.A, [{n}])

        if rec.failed:
            return rec

        # 步骤 (3)-(4)：清理 MD 至 proper
        cur = _cleanup_md(cur, rec, choose)
        if rec.failed:
            return rec

        if check:
            if (v := violates_monotonicity(cur)) is not None:
                rec.failed = f"monotonicity violated: {v}"
                return rec
            cur.check_invariants()

    # 步骤 (6)：MU 已空，逐个切 MD 中的 deg 3 atom
    guard = 0
    while cur.atoms:
        guard += 1
        if guard > max_steps:
            rec.failed = "step limit exceeded in step (6)"
            return rec
        cur = _cleanup_md(cur, rec, choose)
        if rec.failed or not cur.atoms:
            break
        deg3 = [a for a in cur.layer_atoms(Layer.DOWN) if cur.degree(a) == 3]
        if not deg3:
            if not is_proper(cur):
                rec.failed = "MD not proper and no deg 3 atom in step (6)"
                return rec
            remaining = sorted(cur.atoms)
            rec.failed = f"step (6): no deg 3 atom left but atoms remain: {remaining}"
            return rec
        sel = choose.pick("step6", deg3)
        cur = _record(rec, cur, {sel}, OpKind.CUT3, [{sel}])

    return rec
