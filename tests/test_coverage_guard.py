"""枚举覆盖的守卫测试（FINDINGS F11 / LESSONS 17）。

背景：`enumerate_toy1*` 与所有基于单父树的随机生成器都漏掉了「某 MD atom
拥有两个 MD 父」的合法构型。生成器之间互相对拍**发现不了**这个漏洞，因为
它们共享同一个盲区。这里改为拿**定义域的性质**直接质询实现。
"""

import pytest

from molecule_cut.algorithm import run_toy1_algorithm
from molecule_cut.builders import build
from molecule_cut.enumerate import (
    enumerate_toy1_full_labeled,
    enumerate_toy1_general,
    rooted_binary_trees,
)
from molecule_cut.molecule import Layer


def count_two_md_parent_atoms(mol) -> int:
    down = set(mol.layer_atoms(Layer.DOWN))
    return sum(1 for a in down if sum(1 for p in mol.parents(a) if p in down) == 2)


class TestRestrictedGeneratorLimitations:
    """受限枚举器的范围必须明确，不能再被称作完整定义域。"""

    def test_rooted_binary_trees_gives_each_node_one_parent(self):
        """现状：parent 数组每节点恰 1 个父（根为 -1）——这就是盲区来源。"""
        for m in range(1, 7):
            for par in rooted_binary_trees(m):
                assert par[0] == -1
                assert all(0 <= p < i for i, p in enumerate(par) if i > 0)

    @pytest.mark.parametrize("n,nd", [(3, 3), (4, 4), (4, 5)])
    def test_enumerator_currently_omits_two_parent_configs(self, n, nd):
        """现状：枚举器产出的实例中双父 atom 数恒为 0。

            这是受限生成器的设计边界；完整审计改用 `enumerate_toy1_full_labeled`。
        """
        assert all(count_two_md_parent_atoms(m) == 0 for m in enumerate_toy1_general(n, nd))


class TestOmittedConfigsAreLegal:
    """被漏掉的构型确实满足 SPEC，因此漏掉它们是真实的覆盖缺口。"""

    @staticmethod
    def two_parent_instance():
        # MD: d0 -> d1 <- d2（d1 有两个 MD 父）；d0、d2 各挂 2 条跨层 bond
        return build(
            4, 3, [(0, 1), (1, 2), (2, 3)], [(0, 1), (2, 1)], [(0, 0), (1, 0), (2, 2), (3, 2)]
        )

    @staticmethod
    def two_parent_up_instance():
        # MU: u0 -> u2 <- u1；MD 也取一个双父树。两层都满足树和槽位约束。
        return build(
            3, 3, [(0, 2), (1, 2)], [(0, 1), (2, 1)], [(0, 0), (1, 0), (2, 2)]
        )

    def test_instance_satisfies_structural_invariants(self):
        mol = self.two_parent_instance()
        mol.check_invariants()
        assert count_two_md_parent_atoms(mol) == 1

    def test_instance_is_valid_toy_model_one(self):
        """Def 1.2：每个 MU atom 恰有一条跨层 bond。"""
        mol = self.two_parent_instance()
        down = set(mol.layer_atoms(Layer.DOWN))
        for u in mol.layer_atoms(Layer.UP):
            assert sum(1 for c in mol.children(u) if c in down) == 1

    def test_md_is_a_tree(self):
        """Def 1.1(3)：MD 是树——无向意义下连通且边数 = 点数 − 1。"""
        mol = self.two_parent_instance()
        down = set(mol.layer_atoms(Layer.DOWN))
        md_edges = [b for b in mol.bonds if b.parent in down and b.child in down]
        assert len(md_edges) == len(down) - 1
        assert len(mol.components(within=down)) == 1

    def test_algorithm_runs_normally_on_it(self):
        """算法在该实例上正常终止并满足 Prop 1.6。"""
        mol = self.two_parent_instance()
        rec = run_toy1_algorithm(mol)
        assert rec.failed is None, rec.failed
        assert rec.n4 == 1
        assert rec.n33 >= (4 - 1) / 5

    def test_mu_two_parent_instance_is_valid_toy_model_one(self):
        mol = self.two_parent_up_instance()
        mol.check_invariants()
        up = set(mol.layer_atoms(Layer.UP))
        assert any(sum(parent in up for parent in mol.parents(atom)) == 2 for atom in up)
        down = set(mol.layer_atoms(Layer.DOWN))
        assert all(sum(child in down for child in mol.children(atom)) == 1 for atom in up)


def test_full_labeled_enumerator_covers_two_parent_cases_in_both_layers():
    """全域审计器须直接命中曾被两次遗漏的边界构型。"""
    molecules = list(enumerate_toy1_full_labeled(3, 3))
    assert molecules
    up_case = down_case = False
    for mol in molecules:
        up = set(mol.layer_atoms(Layer.UP))
        down = set(mol.layer_atoms(Layer.DOWN))
        up_case |= any(sum(parent in up for parent in mol.parents(atom)) == 2 for atom in up)
        down_case |= any(sum(parent in down for parent in mol.parents(atom)) == 2 for atom in down)
    assert up_case, "完整枚举器未覆盖 MU 双父构型"
    assert down_case, "完整枚举器未覆盖 MD 双父构型"
