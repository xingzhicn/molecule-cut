"""Regression guards for deterministic committed calibration tables."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
GENERATOR_PATH = PROJECT / "scripts" / "generate_data.py"
DATA = PROJECT / "data"


def _load_generator():
    spec = importlib.util.spec_from_file_location("molecule_cut_generate_data", GENERATOR_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _snapshot(directory: Path) -> dict[str, bytes]:
    return {path.name: path.read_bytes() for path in sorted(directory.glob("*.csv"))}


def test_small_generation_is_byte_stable_and_lf_only(tmp_path, monkeypatch):
    generator = _load_generator()
    monkeypatch.setattr(generator, "DATA", tmp_path)

    generator.gen_a_n(3)
    generator.gen_v_d(3, 3)
    generator.gen_coefficients(3)
    first = _snapshot(tmp_path)

    generator.gen_a_n(3)
    generator.gen_v_d(3, 3)
    generator.gen_coefficients(3)
    second = _snapshot(tmp_path)

    assert first == second
    for content in second.values():
        assert b"\r\n" not in content
        header = content.splitlines()[0].decode("utf-8").split(",")
        assert "seconds" not in header


def test_committed_tables_use_canonical_coverage_and_stable_schema():
    expected_rows = {
        "a_n.csv": list(range(1, 8)),
        "v_d.csv": list(range(1, 7)),
    }
    for filename, expected_n in expected_rows.items():
        path = DATA / filename
        content = path.read_bytes()
        assert b"\r\n" not in content
        rows = list(csv.DictReader(content.decode("utf-8").splitlines()))
        assert [int(row["n"]) for row in rows] == expected_n
        assert "seconds" not in rows[0]

    coefficient_text = (DATA / "bond_breaking_coefficients.csv").read_text(encoding="utf-8")
    assert "(6,7)" in coefficient_text
    assert "seconds" not in coefficient_text.splitlines()[0]
