# Reproducibility

This repository is standalone. It uses one Python 3.12 environment managed by
`uv`; no project-local virtual environment is required.

## Verification

```sh
uv sync --frozen
uv run --frozen pytest
uv run --frozen ruff check .
```

The bounded proof-regression audit is:

```sh
uv run --frozen pytest -q tests/test_p_exact_route.py
```

On the current release candidate this focused audit has 4 passing tests; the
full suite has 653 passing tests.  These counts are regression metadata and do
not replace the written all-size proof.

The intrinsic-transfer checks are:

```sh
uv run --frozen pytest -q tests/test_intrinsic.py tests/test_fast_vd.py
```

## Replaying tables

```sh
uv run --frozen python scripts/generate_upper_bounds.py
uv run --frozen python scripts/generate_data.py
```

The default data scripts are deliberately bounded. They are regression and
replay tools, not substitutes for the written all-size proof.

## Finite audit scope

`tests/test_p_exact_route.py` enumerates all labelled Toy I molecules with
positive layer sizes and total size at most six, including slot-legal
orientations, capacity-legal cross-bond assignments, and all tie-break
branches. The frozen audit totals are 11,590 labelled molecules and 25,243
complete tie-break records, with no failed records or branch-cap overflows.
This is finite implementation evidence only; it is not extrapolated to
unrestricted lower-layer size.

## Source

The public source model is Deng, Hani, and Ma, *Long time derivation of the
Boltzmann equation from hard sphere dynamics*, arXiv:2408.07818v3, Sections
11.1--11.2 and Proposition 11.7:
<https://arxiv.org/abs/2408.07818>.
