"""Elementary molecule 分类、proper 判定、monotonicity 检查。

实现 SPEC.md §3.1（elementary 类型）、§3.4（monotonicity）、§3.5（proper）。
"""

from __future__ import annotations

from enum import Enum

from .molecule import Direction, Layer, Molecule


class ElemType(Enum):
    """Overview §3.3 的四种 elementary molecule，外加 {343} 与非法。

    {343} 本身不是 elementary，但 Def 1.4 步骤 (4) 把它切成 {3} + {33}，
    计分时算作一个 {33}（见 (1.3)）。
    """

    E2 = "{2}"  # 单 atom，deg 2，两 free end 同为 top 或同为 bottom
    E3 = "{3}"  # 单 atom，deg 3
    E4 = "{4}"  # 单 atom，deg 4，无 fixed end —— bad
    E33 = "{33}"  # 两 atom 由 bond 相连，各 deg 3 —— good
    E343 = "{343}"  # 两个 deg 3 邻接同一 deg 4；非 elementary，可再切
    INVALID = "invalid"


# 计分权重（SPEC §1.2）：score = #{33} - 10d * #{4}
GOOD_TYPES = frozenset({ElemType.E33})
BAD_TYPES = frozenset({ElemType.E4})


def classify(mol: Molecule, subset: set[int]) -> ElemType:
    """把 subset 作为独立子 molecule 分类（SPEC §1.1 的子 molecule 约定）。"""
    sub = mol.induced(subset)
    n = len(subset)

    if n == 1:
        (aid,) = subset
        deg = sub.degree(aid)
        atom = sub.atoms[aid]
        if deg == 4:
            # {4}: 四个 free end，无 fixed end
            return ElemType.E4 if atom.n_fixed() == 0 else ElemType.INVALID
        if deg == 3:
            return ElemType.E3
        if deg == 2:
            # {2}: 两 free end 必须同为 top 或同为 bottom（Overview §3.3(1)）
            free_dirs = {d for k, d in atom.ends if k.value == "free"}
            return ElemType.E2 if len(free_dirs) == 1 else ElemType.INVALID
        return ElemType.INVALID  # deg 0/1 均非 elementary

    if n == 2:
        a, b = sorted(subset)
        if len(sub.bonds) != 1:
            return ElemType.INVALID
        if sub.degree(a) == 3 and sub.degree(b) == 3:
            return ElemType.E33
        return ElemType.INVALID

    if n == 3:
        # {343}: 中间 deg 4，两端 deg 3，路径结构
        degs = sorted(sub.degree(a) for a in subset)
        if degs != [3, 3, 4]:
            return ElemType.INVALID
        center = [a for a in subset if sub.degree(a) == 4]
        if len(center) != 1:
            return ElemType.INVALID
        c = center[0]
        others = subset - {c}
        # 两个 deg 3 都必须邻接 center，且彼此不相邻（路径而非三角）
        if not all(o in set(sub.neighbors(c)) for o in others):
            return ElemType.INVALID
        o1, o2 = sorted(others)
        if o2 in set(sub.neighbors(o1)):
            return ElemType.INVALID
        return ElemType.E343

    return ElemType.INVALID


def is_elementary(mol: Molecule, subset: set[int]) -> bool:
    return classify(mol, subset) in {
        ElemType.E2,
        ElemType.E3,
        ElemType.E4,
        ElemType.E33,
    }


def is_proper(mol: Molecule, layer: Layer = Layer.DOWN) -> bool:
    """SPEC §3.5 / Def 1.3：MD 为 proper 的三个条件。

    1. 无 deg 1 或 deg 2 atom
    2. 无两个相邻的 deg 3 atom
    3. 无两个 deg 3 atom 同时邻接同一个 deg 4 atom
    """
    ids = mol.layer_atoms(layer)
    deg = {a: mol.degree(a) for a in ids}
    idset = set(ids)

    if any(deg[a] in (1, 2) for a in ids):
        return False

    for a in ids:
        if deg[a] != 3:
            continue
        for nb in mol.neighbors(a):
            if nb in idset and deg[nb] == 3:
                return False  # 两个相邻 deg 3

    for a in ids:
        if deg[a] != 4:
            continue
        deg3_nbrs = [nb for nb in mol.neighbors(a) if nb in idset and deg[nb] == 3]
        if len(deg3_nbrs) >= 2:
            return False  # 两个 deg 3 共享同一 deg 4

    return True


def find_adjacent_deg3_pairs(mol: Molecule, layer: Layer = Layer.DOWN) -> list[set[int]]:
    """Def 1.4 步骤 (3) 的候选：MD 中所有相邻 deg 3 对。tie-breaking 自由度之一。"""
    ids = set(mol.layer_atoms(layer))
    deg = {a: mol.degree(a) for a in ids}
    out: list[set[int]] = []
    for b in mol.bonds:
        p, c = b.parent, b.child
        if p in ids and c in ids and deg[p] == 3 and deg[c] == 3:
            out.append({p, c})
    return out


def find_343_triples(mol: Molecule, layer: Layer = Layer.DOWN) -> list[set[int]]:
    """Def 1.4 步骤 (4) 的候选：两个 deg 3 邻接同一 deg 4。tie-breaking 自由度之一。"""
    ids = set(mol.layer_atoms(layer))
    deg = {a: mol.degree(a) for a in ids}
    out: list[set[int]] = []
    for a in sorted(ids):
        if deg[a] != 4:
            continue
        deg3 = sorted(nb for nb in mol.neighbors(a) if nb in ids and deg[nb] == 3)
        for i in range(len(deg3)):
            for j in range(i + 1, len(deg3)):
                out.append({deg3[i], a, deg3[j]})
    return out


def violates_monotonicity(mol: Molecule) -> str | None:
    """SPEC §3.4：检查是否出现了不可能被切成 elementary 的 atom。

    唯一的非 elementary 单 atom molecule 是 deg 1，或 deg 2 且同时带
    top 和 bottom fixed end（Overview §4.3.1(1)）。返回错误描述或 None。
    """
    for aid, atom in mol.atoms.items():
        deg = mol.degree(aid)
        if deg == 1:
            return f"atom {aid} has deg 1 (forbidden)"
        if deg == 2:
            fd = atom.fixed_directions()
            if Direction.TOP in fd and Direction.BOTTOM in fd:
                return f"atom {aid} has deg 2 with both top and bottom fixed ends"
    return None


def lambda_potential(mol: Molecule, layer: Layer = Layer.DOWN) -> int:
    """SPEC §3.7：λ = #deg3(MD) - #comp(MD)。初值 -1，终值 0。"""
    ids = set(mol.layer_atoms(layer))
    if not ids:
        return 0
    n_deg3 = sum(1 for a in ids if mol.degree(a) == 3)
    n_comp = len(mol.components(within=ids))
    return n_deg3 - n_comp
