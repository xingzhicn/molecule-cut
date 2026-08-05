"""序列化往返测试（支撑 PLAN V1.4：极值实例可反序列化重放）。"""

import pytest

from molecule_cut.algorithm import run_toy1_algorithm
from molecule_cut.enumerate import enumerate_toy1, enumerate_toy1_general
from molecule_cut.exhaustive import a_of_molecule, v_d
from molecule_cut.molecule import Layer
from molecule_cut.serialize import from_dict, from_json, to_dict, to_json


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_roundtrip_preserves_structure(n):
    for mol in enumerate_toy1(n):
        back = from_json(to_json(mol))
        assert back.layer_atoms(Layer.UP) == mol.layer_atoms(Layer.UP)
        assert back.layer_atoms(Layer.DOWN) == mol.layer_atoms(Layer.DOWN)
        assert sorted((b.parent, b.child) for b in back.bonds) == sorted(
            (b.parent, b.child) for b in mol.bonds
        )
        assert back.circuit_rank() == mol.circuit_rank()


@pytest.mark.parametrize("n", [2, 3, 4])
def test_roundtrip_preserves_computed_values(n):
    """重放后 a_n 与 v_d 必须一致——这是 witness 可信的前提。"""
    for mol in enumerate_toy1(n):
        back = from_json(to_json(mol))
        assert a_of_molecule(back).n33_min == a_of_molecule(mol).n33_min
        assert v_d(back, d=3) == v_d(mol, d=3)
        assert run_toy1_algorithm(back).n33 == run_toy1_algorithm(mol).n33


@pytest.mark.parametrize("n,nd", [(3, 2), (3, 4), (4, 3), (4, 5)])
def test_roundtrip_with_unequal_layers(n, nd):
    for mol in enumerate_toy1_general(n, nd):
        back = from_dict(to_dict(mol))
        assert len(back.atoms) == len(mol.atoms)
        assert len(back.bonds) == len(mol.bonds)


def test_to_dict_rejects_noncanonical_ids():
    """切割后的 molecule 编号不连续，应明确报错而不是静默产出错误数据。"""
    mol = next(enumerate_toy1(3))
    cut = mol.cut_as_free({0})
    with pytest.raises(ValueError, match="canonical ids"):
        to_dict(cut)


def test_from_dict_revalidates_invariants():
    """反序列化必须重新执行构造期检查（槽位越界应被拒绝）。"""
    bad = {
        "up": 1,
        "down": 1,
        "up_edges": [],
        "down_edges": [],
        "cross": [[0, 0], [0, 0], [0, 0]],  # 同一对之间 3 条 bond，超出槽位
    }
    with pytest.raises(ValueError):
        from_dict(bad)
