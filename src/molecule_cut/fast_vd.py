"""v_d(M) 的快速实现（位掩码状态 + 预计算静态结构）。

**关键观察**：切割过程中的状态**完全由「剩余 atom 集合」决定**。
因为 atom a 的 fixed end 数 = 其被移除的原邻居数，free end 数不变，故

    deg(a, S) = |nbrs(a) ∩ S| + base_free(a),   base_free(a) = 4 − |nbrs(a)|

于是状态空间至多 2^|M|，而 `exhaustive.v_d` 的瓶颈不是状态数，而是每状态开销：
它为每个候选切割构造一个新的 `Molecule`（`induced()`），并反复访问 Enum 属性。

本模块把状态表示为整数位掩码，静态结构（邻居、父子、槽位）一次性预计算，
分类不构造任何对象。语义与 `exhaustive.v_d` 完全一致，由差分测试保证
（`tests/test_fast_vd.py`）。
"""

from __future__ import annotations

from .molecule import Layer, Molecule

# 组件类型的轻量编码（避免热路径上的 Enum 访问）
BAD4 = 0  # {4}
NORM3 = 1  # {3}
NORM2 = 2  # {2}
GOOD33 = 3  # {33}
GOOD343 = 4  # {343}，拆成 {3} + {33}，计一个好分量
INVALID = -1


class FastMolecule:
    """预计算的静态结构，供位掩码状态查询。"""

    __slots__ = (
        "base_free",
        "chi_mask",
        "full_mask",
        "ids",
        "index",
        "is_up",
        "n_atoms",
        "n_orig_chi",
        "n_orig_nbr",
        "n_orig_par",
        "nbr_mask",
        "par_mask",
    )

    def __init__(self, mol: Molecule) -> None:
        self.ids = sorted(mol.atoms)
        self.index = {a: i for i, a in enumerate(self.ids)}
        n = self.n_atoms = len(self.ids)
        self.full_mask = (1 << n) - 1

        self.nbr_mask = [0] * n
        self.par_mask = [0] * n
        self.chi_mask = [0] * n
        self.n_orig_par = [0] * n
        self.n_orig_chi = [0] * n
        self.is_up = [mol.atoms[a].layer is Layer.UP for a in self.ids]

        for b in mol.bonds:
            p, c = self.index[b.parent], self.index[b.child]
            self.nbr_mask[p] |= 1 << c
            self.nbr_mask[c] |= 1 << p
            self.chi_mask[p] |= 1 << c
            self.par_mask[c] |= 1 << p
            self.n_orig_chi[p] += 1
            self.n_orig_par[c] += 1

        self.n_orig_nbr = [self.n_orig_par[i] + self.n_orig_chi[i] for i in range(n)]
        # 每个 atom 共 4 个槽位；未被 bond 占用的即初始 free end
        self.base_free = [4 - self.n_orig_nbr[i] for i in range(n)]

    def degree(self, i: int, state: int) -> int:
        return (self.nbr_mask[i] & state).bit_count() + self.base_free[i]

    def n_fixed(self, i: int, state: int) -> int:
        """被移除的原邻居数。"""
        return self.n_orig_nbr[i] - (self.nbr_mask[i] & state).bit_count()

    def classify(self, subset: list[int], state: int) -> int:
        """给 subset（状态 state 下的下标列表）分类，语义同 elementary.classify。"""
        k = len(subset)
        if k == 1:
            i = subset[0]
            deg = self.degree(i, state)
            if deg == 4:
                return BAD4 if self.n_fixed(i, state) == 0 else INVALID
            if deg == 3:
                return NORM3
            if deg == 2:
                # 切出后 top/bottom 两侧的 free end 数；{2} 要求二者只有一侧非零
                top = (2 - self.n_orig_par[i]) + (self.par_mask[i] & state).bit_count()
                bot = (2 - self.n_orig_chi[i]) + (self.chi_mask[i] & state).bit_count()
                return NORM2 if (top == 0) != (bot == 0) else INVALID
            return INVALID

        if k == 2:
            a, b = subset
            if not (self.nbr_mask[a] >> b) & 1:
                return INVALID
            if self.degree(a, state) == 3 and self.degree(b, state) == 3:
                return GOOD33
            return INVALID

        if k == 3:
            degs = [self.degree(i, state) for i in subset]
            if sorted(degs) != [3, 3, 4]:
                return INVALID
            centre = subset[degs.index(4)]
            others = [i for i in subset if i != centre]
            for o in others:
                if not (self.nbr_mask[centre] >> o) & 1:
                    return INVALID
            if (self.nbr_mask[others[0]] >> others[1]) & 1:
                return INVALID  # 三角形而非路径
            return GOOD343

        return INVALID

    def legal_moves(self, state: int) -> list[tuple[tuple[int, ...], int]]:
        """当前状态下所有合法切割及其类型。"""
        out: list[tuple[tuple[int, ...], int]] = []
        alive = [i for i in range(self.n_atoms) if (state >> i) & 1]

        for i in alive:
            t = self.classify([i], state)
            if t != INVALID:
                out.append(((i,), t))

        for i in alive:
            nb = self.nbr_mask[i] & state
            j = i + 1
            m = nb >> j
            while m:
                if m & 1 and self.classify([i, j], state) == GOOD33:
                    out.append(((i, j), GOOD33))
                m >>= 1
                j += 1

        for c in alive:
            if self.degree(c, state) != 4:
                continue
            deg3 = [
                i
                for i in alive
                if (self.nbr_mask[c] >> i) & 1 and self.degree(i, state) == 3
            ]
            for x in range(len(deg3)):
                for y in range(x + 1, len(deg3)):
                    trip = [deg3[x], c, deg3[y]]
                    if self.classify(trip, state) == GOOD343:
                        out.append((tuple(sorted(trip)), GOOD343))

        return out


def _solve(fm: FastMolecule, gain_of: dict[int, int]) -> int | None:
    """在固定的逐步得分下做位掩码 DP。"""
    memo: dict[int, int | None] = {}

    def rec(state: int) -> int | None:
        if state == 0:
            return 0
        cached = memo.get(state, ...)
        if cached is not ...:
            return cached  # type: ignore[return-value]
        best: int | None = None
        for subset, t in fm.legal_moves(state):
            nxt = state
            for i in subset:
                nxt &= ~(1 << i)
            sub = rec(nxt)
            if sub is None:
                continue
            val = gain_of[t] + sub
            if best is None or val > best:
                best = val
        memo[state] = best
        return best

    return rec(fm.full_mask)


def fast_v_d(mol: Molecule, d: int = 3) -> int | None:
    """计算 `v_d`，即 `#{33} - 10d·#{4}` 的最优值。"""
    fm = FastMolecule(mol)
    return _solve(fm, {BAD4: -10 * d, NORM3: 0, NORM2: 0, GOOD33: 1, GOOD343: 1})


def fast_w(mol: Molecule) -> int | None:
    """计算 `W(M) = max_C #{33}(C)`，不把 `{4}` 惩罚混入目标。"""
    fm = FastMolecule(mol)
    return _solve(fm, {BAD4: 0, NORM3: 0, NORM2: 0, GOOD33: 1, GOOD343: 1})
