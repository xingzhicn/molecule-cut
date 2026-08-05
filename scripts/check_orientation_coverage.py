#!/usr/bin/env python3
"""检查枚举覆盖漏洞：MD 中允许 atom 拥有两个 MD 父。

背景（2026-08-01 与 2026-08-04 审查发现）：受限的 `enumerate_toy1*` 固定 MU
为链、MD 为单父有根树；但 Def 1.1(3) 只要求两层是**树**。定向可使 MU 或 MD
atom 有两个 layer-internal parent，且这类构型合法（`build` 接受、
`check_invariants` 通过、算法可正常运行）。

本脚本用完整的「两层无向树 × 所有合法定向」标号枚举重算 a_n，与受限枚举对照。
它只适合很小规模；不再把任何 n=5 以上输出称为已经完成的全枚举。

用法：
    uv run --frozen python scripts/check_orientation_coverage.py --max-n 3
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from molecule_cut.algorithm import run_toy1_algorithm
from molecule_cut.enumerate import enumerate_toy1, enumerate_toy1_full_labeled
from molecule_cut.exhaustive import a_of_molecule
from molecule_cut.molecule import Layer
from molecule_cut.serialize import to_json


def enumerate_full(n: int, md: int):
    """|MU| = n、|MD| = md 的完整标号 Toy I 枚举。

    不做同构去重（只用于与受限枚举比较最小值，去重不影响 min）。
    """
    yield from enumerate_toy1_full_labeled(n, md)


def has_two_md_parents(mol) -> bool:
    down = set(mol.layer_atoms(Layer.DOWN))
    return any(sum(1 for p in mol.parents(a) if p in down) == 2 for a in down)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-n", type=int, default=3)
    args = p.parse_args()

    print("n  受限枚举 a_n  完整枚举 a_n  ⌈(n-1)/3⌉  含双父实例数  结论")
    for n in range(2, args.max_n + 1):
        restricted = min(a_of_molecule(m).n33_min for m in enumerate_toy1(n))

        full_best = None
        witness = None
        two_parent_count = 0
        for md in range(max(1, n - 1), n + 2):
            for mol in enumerate_full(n, md):
                if has_two_md_parents(mol):
                    two_parent_count += 1
                rec = run_toy1_algorithm(mol, check=False)
                if rec.failed:
                    continue
                v = a_of_molecule(mol).n33_min
                if full_best is None or v < full_best:
                    full_best, witness = v, mol

        target = math.ceil((n - 1) / 3)
        verdict = "一致" if full_best == restricted else f"**不同！完整枚举更小 {full_best}**"
        if full_best is not None and full_best < target:
            verdict += "  且低于 F7 预测 —— 反例"
        print(
            f"{n}  {restricted:11d}  {full_best!s:>11s}  {target:8d}  {two_parent_count:11d}  {verdict}"
        )
        if full_best is not None and full_best < restricted and witness is not None:
            print(f"    witness: {to_json(witness)}")


if __name__ == "__main__":
    main()
