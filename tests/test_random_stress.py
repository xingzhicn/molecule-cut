"""随机 toy model I 上的不变量压力测试。

固定种子以保证可复现（SPEC 无随机性要求，这里只是扩大覆盖面）。
"""

import random

import pytest

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.builders import build

SEED = 20260728
N_SAMPLES = 300


def random_binary_forest_edges(n: int, rng: random.Random) -> list[tuple[int, int]]:
    """随机生成 n 个节点的树（每节点至多 2 child），返回 (parent, child) 列表。"""
    edges: list[tuple[int, int]] = []
    child_count = {0: 0}
    for v in range(1, n):
        avail = [u for u in range(v) if child_count[u] < 2]
        p = rng.choice(avail)
        edges.append((p, v))
        child_count[p] += 1
        child_count[v] = 0
    return edges


def random_toy1(n: int, rng: random.Random):
    """随机 toy model I：MU 为链（保证唯一最低 atom），MD 为随机二叉树。

    MU 用链只是最简单的合法抽样方式；Def 1.2 也允许有双父 atom 的倒置二叉树。
    """
    up_edges = [(i, i + 1) for i in range(n - 1)]
    down_edges = random_binary_forest_edges(n, rng)
    cross = [(i, p) for i, p in enumerate(rng.sample(range(n), n))]
    return build(n, n, up_edges, down_edges, cross)


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6, 7, 8])
def test_random_toy1_invariants(n):
    rng = random.Random(SEED + n)
    for trial in range(N_SAMPLES // n):
        try:
            mol = random_toy1(n, rng)
        except ValueError:
            continue  # 槽位冲突（某 atom 超过 2 parent），跳过
        rec = run_toy1_algorithm(mol)
        ctx = f"n={n} trial={trial}"

        assert rec.failed is None, f"{ctx}: {rec.failed}"

        c = rec.op_counts()
        n_a = c.get(OpKind.A, 0)
        n_b = c.get(OpKind.B, 0)
        n_33 = c.get(OpKind.CUT33, 0)
        n_343 = c.get(OpKind.CUT343, 0)

        # (1.3)
        assert rec.n33 == n_b + n_33 + n_343, f"{ctx}: identity (1.3)"
        # (1.4)
        assert n_a - (n_33 + n_343) == 1, f"{ctx}: identity (1.4)"
        # bond counting
        assert n_a + 2 * n_b + 3 * n_33 + 4 * n_343 >= n, f"{ctx}: bond counting"
        # Prop 1.6
        assert rec.n4 == 1, f"{ctx}: #{{4}} != 1"
        assert rec.n33 >= (n - 1) / 5, f"{ctx}: #{{33}} below (n-1)/5"
        # λ 端点
        assert rec.steps[0].lam_before == -1, f"{ctx}: lambda start"
        assert rec.steps[-1].lam_after == 0, f"{ctx}: lambda end"
        # atom 覆盖
        cut = [a for s in rec.steps for a in s.subset]
        assert sorted(cut) == sorted(mol.atoms), f"{ctx}: atom coverage"
