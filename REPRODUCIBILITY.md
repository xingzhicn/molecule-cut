# Reproducibility

This is a standalone Python 3.12 project managed by [uv](https://docs.astral.sh/uv/).
The committed `uv.lock` is part of the release contract.

## Verification

From the repository root, run:

```sh
uv sync --frozen --group dev
uv run --frozen pytest
uv run --frozen ruff check .
```

GitHub Actions executes the same test and lint commands on pushes to `main` and
on pull requests. The tests are regression and transcription checks; the
all-size theorem is proved in `MANUSCRIPT.pdf`, not extrapolated from them.

For release `v0.1.0`, these commands reported 675 passing tests and a clean
Ruff check.

The focused finite proof-regression audit is:

```sh
uv run --frozen pytest -q tests/test_p_exact_route.py
```

For release `v0.1.0`, this focused audit reported 6 passing tests.

The intrinsic-objective checks are:

```sh
uv run --frozen pytest -q tests/test_intrinsic.py tests/test_fast_vd.py
```

## Replaying data

```sh
uv run --frozen python scripts/generate_upper_bounds.py
uv run --frozen python scripts/generate_data.py
```

The stored tables are deterministic bounded computations. `generate_upper_bounds.py`
replays the explicit sharpness families; `generate_data.py` replays the broader
calibration tables. Neither command is a premise of the infinite theorem.

## Finite audit scope

`tests/test_p_exact_route.py` enumerates every labelled Toy I molecule with
positive layer sizes and total size at most six, including slot-legal
orientations, capacity-legal cross-bond assignments, and every tie-break branch.
Its frozen totals are 11,590 labelled molecules and 25,243 complete tie-break
records, with no failed records or branch-cap overflows. The same module gives
complete labelled `(3,3)` and `(4,3)` counterexamples to global tie-invariance.

These are finite implementation checks only. They are not extrapolated to
unrestricted lower-layer size or used as a premise of the proof.

## Source

The public source model is Deng, Hani, and Ma, *Long time derivation of the
Boltzmann equation from hard sphere dynamics*, arXiv:2408.07818v3, Sections
11.1--11.2 and Proposition 11.7:
<https://arxiv.org/abs/2408.07818>.
