"""Molecule 数据结构与 cut 操作。

严格实现 SPEC.md §1（对象定义）与 §3（cut 语义）。术语为纯组合术语，
不含 PDE 内容（讲义 Overview §1.1 的明确要求）。

关键约定（SPEC §3.1）：degree = bond 数 + free end 数，**不计 fixed end**。
每个 atom 的 bonds + free + fixed 恒为 4。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

ATOM_SLOTS = 4  # 每个 atom 的总端点数（bond + free + fixed），讲义约定


class Layer(Enum):
    UP = "U"
    DOWN = "D"


class EndKind(Enum):
    """半边的三种状态。bond 是与另一 atom 相连；free/fixed 是悬空端。"""

    FREE = "free"
    FIXED = "fixed"


class Direction(Enum):
    """半边的方向。bond (p, q) 中 p 是 q 的 parent 时，该边在 q 处是 TOP 方向。"""

    TOP = "top"
    BOTTOM = "bottom"


@dataclass(frozen=True)
class Bond:
    """一条 bond。parent 在上，child 在下。"""

    parent: int
    child: int

    def other(self, atom: int) -> int:
        if atom == self.parent:
            return self.child
        if atom == self.child:
            return self.parent
        raise ValueError(f"atom {atom} not in bond {self}")

    def normalized(self) -> tuple[int, int]:
        return (self.parent, self.child)


@dataclass
class Atom:
    """一个 atom。

    ends 记录非 bond 的半边：(kind, direction) 的多重集合，用列表表示。
    bond 由 Molecule 统一维护，不在 Atom 内冗余存储。
    """

    aid: int
    layer: Layer
    ends: list[tuple[EndKind, Direction]] = field(default_factory=list)

    def n_free(self) -> int:
        return sum(1 for k, _ in self.ends if k is EndKind.FREE)

    def n_fixed(self) -> int:
        return sum(1 for k, _ in self.ends if k is EndKind.FIXED)

    def fixed_directions(self) -> set[Direction]:
        return {d for k, d in self.ends if k is EndKind.FIXED}


class Molecule:
    """两层 molecule。

    不变量（构造后与每次 cut 后均应成立，见 check_invariants）：
      * 每个 atom：len(bonds_at(a)) + n_free + n_fixed == 4
      * 每个 atom：至多 2 个 parent、至多 2 个 child
    """

    def __init__(self, atoms: dict[int, Atom], bonds: list[Bond]) -> None:
        self.atoms = atoms
        self.bonds = list(bonds)

    # ---------- 基本查询 ----------

    def bonds_at(self, aid: int) -> list[Bond]:
        return [b for b in self.bonds if b.parent == aid or b.child == aid]

    def degree(self, aid: int) -> int:
        """SPEC §3.1：deg = bond 数 + free end 数，fixed end 不计。"""
        atom = self.atoms[aid]
        return len(self.bonds_at(aid)) + atom.n_free()

    def parents(self, aid: int) -> list[int]:
        return [b.parent for b in self.bonds if b.child == aid]

    def children(self, aid: int) -> list[int]:
        return [b.child for b in self.bonds if b.parent == aid]

    def layer_atoms(self, layer: Layer) -> list[int]:
        return sorted(a for a, at in self.atoms.items() if at.layer is layer)

    def cross_bonds(self) -> list[Bond]:
        """跨层 bond（parent 在 UP、child 在 DOWN；Def 1.1(2) 保证方向唯一）。"""
        return [
            b for b in self.bonds if self.atoms[b.parent].layer is not self.atoms[b.child].layer
        ]

    def neighbors(self, aid: int) -> list[int]:
        return [b.other(aid) for b in self.bonds_at(aid)]

    # ---------- 连通性 ----------

    def components(self, within: set[int] | None = None) -> list[set[int]]:
        """无向连通分支。within 限定考察的 atom 子集（默认全部）。"""
        pool = set(self.atoms) if within is None else set(within)
        seen: set[int] = set()
        comps: list[set[int]] = []
        for start in sorted(pool):
            if start in seen:
                continue
            stack, comp = [start], set()
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                for nb in self.neighbors(cur):
                    if nb in pool and nb not in comp:
                        stack.append(nb)
            seen |= comp
            comps.append(comp)
        return comps

    def circuit_rank(self) -> int:
        """SPEC §1.3：独立环数 = |E| - |V| + #components。"""
        n_comp = len(self.components())
        return len(self.bonds) - len(self.atoms) + n_comp

    # ---------- cut 操作 ----------

    def cut_as_free(self, subset: set[int]) -> Molecule:
        """SPEC §3.2：切 subset 为 free，返回新 Molecule（不修改自身）。

        对每条跨界 bond (p in subset, q outside)：
          * q 处该 bond 变为 fixed end，方向按 p 相对 q 的位置确定
          * q 的 deg 因此减 1（fixed 不计度）
        subset 及其内部 bond 一并移除。
        """
        if not subset <= set(self.atoms):
            raise ValueError("subset contains unknown atoms")

        new_atoms = {
            aid: Atom(aid, at.layer, list(at.ends))
            for aid, at in self.atoms.items()
            if aid not in subset
        }
        new_bonds: list[Bond] = []

        for b in self.bonds:
            p_in, c_in = b.parent in subset, b.child in subset
            if p_in and c_in:
                continue  # 内部 bond，随 subset 一并移除
            if not p_in and not c_in:
                new_bonds.append(b)
                continue
            # 跨界：外部端点获得一个 fixed end
            if p_in:
                # parent 被切走 → 在 child 处留下 top fixed end
                new_atoms[b.child].ends.append((EndKind.FIXED, Direction.TOP))
            else:
                # child 被切走 → 在 parent 处留下 bottom fixed end
                new_atoms[b.parent].ends.append((EndKind.FIXED, Direction.BOTTOM))

        return Molecule(new_atoms, new_bonds)

    def induced(self, subset: set[int]) -> Molecule:
        """把 subset 视作独立 molecule（SPEC §1.1 子 molecule 约定：
        与外部的 bond 在其中视为 free end）。用于给切出的组件分类。
        """
        atoms = {
            aid: Atom(aid, self.atoms[aid].layer, list(self.atoms[aid].ends)) for aid in subset
        }
        bonds: list[Bond] = []
        for b in self.bonds:
            p_in, c_in = b.parent in subset, b.child in subset
            if p_in and c_in:
                bonds.append(b)
            elif p_in:
                atoms[b.parent].ends.append((EndKind.FREE, Direction.BOTTOM))
            elif c_in:
                atoms[b.child].ends.append((EndKind.FREE, Direction.TOP))
        return Molecule(atoms, bonds)

    # ---------- 不变量 ----------

    def check_invariants(self) -> None:
        """结构不变量。任何违反都说明实现或输入有误。"""
        for aid, atom in self.atoms.items():
            total = len(self.bonds_at(aid)) + len(atom.ends)
            if total != ATOM_SLOTS:
                raise AssertionError(
                    f"atom {aid}: slots={total} != {ATOM_SLOTS} "
                    f"(bonds={len(self.bonds_at(aid))}, ends={len(atom.ends)})"
                )
            if len(self.parents(aid)) > 2:
                raise AssertionError(f"atom {aid} has >2 parents")
            if len(self.children(aid)) > 2:
                raise AssertionError(f"atom {aid} has >2 children")

    def copy(self) -> Molecule:
        return Molecule(
            {aid: Atom(aid, at.layer, list(at.ends)) for aid, at in self.atoms.items()},
            list(self.bonds),
        )

    def __repr__(self) -> str:
        return (
            f"Molecule(|U|={len(self.layer_atoms(Layer.UP))}, "
            f"|D|={len(self.layer_atoms(Layer.DOWN))}, "
            f"bonds={len(self.bonds)}, rho={self.circuit_rank()})"
        )
