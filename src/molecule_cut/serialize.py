"""Molecule 的序列化与反序列化。

用途：把极值实例随数据表一起保存，使任何结论都能被独立重放验证
（PLAN V1.4：「每个极值附 molecule 实例（可反序列化重放）」）。

格式选用 builders.build 的入参形式（up/down/up_edges/down_edges/cross），
而不是直接 dump 内部状态——这样反序列化会重新走一遍构造期不变量检查。
"""

from __future__ import annotations

import json

from .builders import build
from .molecule import Layer, Molecule


def to_dict(mol: Molecule) -> dict:
    """导出为 builders.build 的入参形式。

    要求 mol 是未经切割的初始 molecule（MU 编号 0..up-1，MD 编号 up..up+down-1）。
    """
    up_ids = mol.layer_atoms(Layer.UP)
    down_ids = mol.layer_atoms(Layer.DOWN)
    up, down = len(up_ids), len(down_ids)

    if up_ids != list(range(up)) or down_ids != list(range(up, up + down)):
        raise ValueError("to_dict expects canonical ids: MU = 0..up-1, MD = up..up+down-1")

    up_edges, down_edges, cross = [], [], []
    for b in mol.bonds:
        p_up = b.parent < up
        c_up = b.child < up
        if p_up and c_up:
            up_edges.append([b.parent, b.child])
        elif not p_up and not c_up:
            down_edges.append([b.parent - up, b.child - up])
        elif p_up and not c_up:
            cross.append([b.parent, b.child - up])
        else:
            raise ValueError(f"bond {b} goes from MD to MU, violating Def 1.1(2)")

    return {
        "up": up,
        "down": down,
        "up_edges": sorted(up_edges),
        "down_edges": sorted(down_edges),
        "cross": sorted(cross),
    }


def from_dict(data: dict) -> Molecule:
    """从 to_dict 的输出重建 molecule（会重新执行构造期不变量检查）。"""
    return build(
        up=data["up"],
        down=data["down"],
        up_edges=[tuple(e) for e in data["up_edges"]],
        down_edges=[tuple(e) for e in data["down_edges"]],
        cross=[tuple(e) for e in data["cross"]],
    )


def to_json(mol: Molecule) -> str:
    """紧凑单行 JSON，便于放进 csv 单元格。"""
    return json.dumps(to_dict(mol), separators=(",", ":"), sort_keys=True)


def from_json(s: str) -> Molecule:
    return from_dict(json.loads(s))
