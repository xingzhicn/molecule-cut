#!/usr/bin/env python3
"""重建枚举、DP 与系数数据表（PLAN V1.6）。

用法（在课题目录下）：

    uv run --frozen python scripts/generate_data.py            # 默认规模
    uv run --frozen python scripts/generate_data.py --max-n 6  # 限制规模

只依赖标准库与本包，不需要 research 依赖组。

本脚本产出：
  data/a_n.csv                        a_n 精确值 + 极值 witness
  data/v_d.csv                        W(M)/v_d 的范围 + min-max 间隙 + witness
  data/bond_breaking_coefficients.csv 每种操作实测能打断的跨层 bond 数

F9 的构造性上界表由独立命令
`uv run --frozen python scripts/generate_upper_bounds.py` 重建。
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from molecule_cut.algorithm import OpKind
from molecule_cut.enumerate import enumerate_toy1, enumerate_toy1_general
from molecule_cut.exhaustive import (
    a_of_molecule,
    enumerate_tiebreak_records,
    v_d,  # 默认 fast path
    w,
)
from molecule_cut.serialize import to_json

DATA = ROOT / "data"
D_DEFAULT = 3

# 讲义 Prop 1.6 证明 part 4 使用的系数（SPEC §3.7）
NOTES_COEFFICIENT = {OpKind.A: 1, OpKind.B: 2, OpKind.CUT33: 3, OpKind.CUT343: 4}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.relative_to(ROOT)} ({len(rows)} rows)")


def gen_a_n(max_n: int) -> None:
    """a_n = min over M, min over tie-breaks, #{33}。

    精确性有两种来源，均给出精确值，在 exactness 列中区分：
      * 找到达到 Prop 1.6 已证下界 ⌈(n-1)/5⌉ 的实例 → 上下界相遇
      * 全枚举
    """
    print("a_n ...")
    rows = []
    for n in range(1, max_n + 1):
        lb = math.ceil((n - 1) / 5)
        t = time.time()
        best = None
        witness = None
        scanned = 0
        for mol in enumerate_toy1(n):
            r = a_of_molecule(mol)
            scanned += 1
            if best is None or r.n33_min < best:
                best, witness = r.n33_min, mol
            if best <= lb:
                break
        hit_lb = best <= lb
        rows.append(
            {
                "n": n,
                "a_n": best,
                "prop16_lower_bound": lb,
                "exactness": (
                    "exact_hit_proven_lower_bound"
                    if hit_lb
                    else "exact_full_enumeration_of_restricted_subclass"
                ),
                "molecules_scanned": scanned,
                "seconds": round(time.time() - t, 2),
                # 此枚举器固定 MU 为链、MD 为单父有根树；不是完整 Toy I 定义域。
                "enumeration_scope": "|MD|=|MU|=n; MU chain; MD single-parent rooted tree",
                "witness_json": to_json(witness),
            }
        )
        print(f"  n={n}: a_n={best} (lb={lb}, scanned={scanned})", flush=True)
    write_csv(DATA / "a_n.csv", rows)


def gen_v_d(max_n: int, d: int) -> None:
    """独立计算 W 与 v_d，并记录二者的 min-max 数据。"""
    print("v_d / W ...")
    rows = []
    for n in range(1, max_n + 1):
        t = time.time()
        ws, vs, ans = [], [], []
        w_min_witness = None
        w_min = None
        for mol in enumerate_toy1(n):
            v = v_d(mol, d=d)
            if v is None:
                raise RuntimeError(f"n={n}: molecule with no legal complete cut: {to_json(mol)}")
            value_w = w(mol)
            if value_w is None:
                raise RuntimeError(f"n={n}: molecule with no legal complete cutting for W: {to_json(mol)}")
            ws.append(value_w)
            vs.append(v)
            ans.append(a_of_molecule(mol).n33_min)
            if w_min is None or value_w < w_min:
                w_min, w_min_witness = value_w, mol
        rows.append(
            {
                "n": n,
                "rho": n - 1,
                "d": d,
                "W_min": min(ws),
                "W_max": max(ws),
                "v_d_min": min(vs),
                "v_d_max": max(vs),
                "a_n": min(ans),
                "minmax_gap": min(ws) - min(ans),
                "molecules": len(ws),
                "seconds": round(time.time() - t, 1),
                "enumeration_scope": "|MD|=|MU|=n; MU chain; MD single-parent rooted tree",
                "w_min_witness_json": to_json(w_min_witness),
            }
        )
        print(f"  n={n}: W_min={min(ws)} a_n={min(ans)} gap={min(ws) - min(ans)}", flush=True)
    write_csv(DATA / "v_d.csv", rows)


def gen_coefficients(max_n: int) -> None:
    """实测每种操作能打断多少条跨层 bond，对比讲义证明使用的系数。

    这是 FINDINGS F3 的数据来源：讲义用 (1,2,3,4)，实测上界更小则可收紧 LP。
    """
    print("bond-breaking coefficients ...")
    mx: dict[OpKind, int] = defaultdict(int)
    hist: dict[OpKind, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    scope = []
    for n in range(3, max_n + 1):
        for nd in (n - 1, n, n + 1):
            if nd < 1:
                continue
            scope.append(f"({n},{nd})")
            for mol in enumerate_toy1_general(n, nd):
                for _, rec in enumerate_tiebreak_records(mol):
                    if rec.failed:
                        continue
                    for s in rec.steps:
                        mx[s.op] = max(mx[s.op], s.cross_bonds_broken)
                        hist[s.op][s.cross_bonds_broken] += 1
    # LESSONS 18 / FINDINGS F3 + F11：这里的 observed_max **不是真实上界**。
    # 本枚举器只覆盖 MU 为链、MD 为单父有根树的子类；在完整定义域上 cut{33} 可达 3、
    # cut{343} 可达 4（反例见 tests/test_coefficient_witnesses.py）。
    scope_str = (
        f"(|MU|,|MD|) in {{{','.join(scope)}}}, all tie-breaks; "
        "MU restricted to a chain; MD restricted to single-parent rooted trees (F11) -- "
        "observed_max is NOT a true upper bound, see F3"
    )
    rows = []
    for op in OpKind:
        if op not in mx:
            continue
        claimed = NOTES_COEFFICIENT.get(op)
        rows.append(
            {
                "operation": op.value,
                "notes_coefficient": claimed if claimed is not None else "",
                "observed_max": mx[op],
                "gap": (claimed - mx[op]) if claimed is not None else "",
                "distribution": str(dict(sorted(hist[op].items()))),
                "scope": scope_str,
            }
        )
        print(f"  {op.value:10s} notes={claimed} observed_max={mx[op]}", flush=True)
    write_csv(DATA / "bond_breaking_coefficients.csv", rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-n", type=int, default=5, help="a_n 的最大 n（默认 5）")
    p.add_argument(
        "--max-n-vd",
        type=int,
        default=5,
        help="v_d/W 的最大 n（默认 5；两种目标均独立计算）",
    )
    p.add_argument("--max-n-coef", type=int, default=3, help="系数实测的最大 n（默认 3）")
    p.add_argument("--d", type=int, default=D_DEFAULT, help="维数参数 d（默认 3）")
    args = p.parse_args()

    t0 = time.time()
    gen_a_n(args.max_n)
    gen_v_d(args.max_n_vd, args.d)
    gen_coefficients(args.max_n_coef)
    print(f"done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
