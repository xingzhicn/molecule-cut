#!/usr/bin/env python3
"""重建 data/a_n_upper_bounds.csv（PLAN V1.6：单条命令可复现）。

F9 的显式族 `family_for(n)` 对每个 n ≥ 4 直接构造达到 ⌈(n−1)/3⌉ 的实例，
故上界是**构造性**的，不再依赖随机搜索。（此前 n = 8..18 的行来自随机搜索，
其中 n = 18 因搜索不足只得到 7 > 6；本脚本一并修正。）

用法（在课题目录下）：

    uv run --frozen python scripts/generate_upper_bounds.py
    uv run --frozen python scripts/generate_upper_bounds.py --max-n 60
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from molecule_cut.algorithm import OpKind, run_toy1_algorithm
from molecule_cut.families import family_for
from molecule_cut.molecule import Layer
from molecule_cut.serialize import to_json

DATA = ROOT / "data"
CLEANUP = (OpKind.CUT33, OpKind.CUT343)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-n", type=int, default=40)
    args = p.parse_args()

    rows = []
    for n in range(4, args.max_n + 1):
        mol = family_for(n)
        mol.check_invariants()
        rec = run_toy1_algorithm(mol)
        if rec.failed:
            raise RuntimeError(f"n={n}: family run failed: {rec.failed}")

        target = math.ceil((n - 1) / 3)
        steps = [s for s in rec.steps if s.op in CLEANUP]
        g = len(steps)
        sum_b = sum(s.cross_bonds_broken for s in steps)

        rows.append(
            {
                "n": n,
                "a_n_upper_bound": rec.n33,
                "ceil_n_minus_1_over_3": target,
                "prop16_lower_bound": math.ceil((n - 1) / 5),
                "attains_target": rec.n33 == target,
                "family": {1: "transpose", 2: "padded", 0: "trimmed"}[n % 3],
                "md_size": len(mol.layer_atoms(Layer.DOWN)),
                "cleanup_ops_g": g,
                "sum_b_i": sum_b,
                "two_g": 2 * g,
                "provenance": "constructive: families.family_for(n), FINDINGS F9",
                "witness_json": to_json(mol),
            }
        )

    bad = [r["n"] for r in rows if not r["attains_target"]]
    if bad:
        raise RuntimeError(f"family failed to attain target at n={bad}")

    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / "a_n_upper_bounds.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out.relative_to(ROOT)} ({len(rows)} rows, n = 4..{args.max_n})")
    print(
        f"  all attain ceil((n-1)/3); sum_b_i <= 2g on every row: "
        f"{all(r['sum_b_i'] <= r['two_g'] for r in rows)}"
    )


if __name__ == "__main__":
    main()
