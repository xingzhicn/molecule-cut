# molecule-cut

[![Verify](https://github.com/xingzhicn/molecule-cut/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/xingzhicn/molecule-cut/actions/workflows/ci.yml)

Reproducible code and preprint for the sharp Toy I benchmark in the
Deng--Hani--Ma two-layer molecule cutting problem.

## Main result

Let $N_{33}$ denote the number of `{33}` components (the source's `#{33}`
notation), and let $n=|M_U|$. For every finite Toy I molecule and every
permitted tie-break of the source Definition 11.4 algorithm, the paper proves

$$
N_{33}\ge \left\lceil\frac{n-1}{3}\right\rceil.
$$

This bound is sharp:

$$
a_n=\left\lceil\frac{n-1}{3}\right\rceil \qquad (n\ge4),
$$

where $a_n$ is the minimum number of `{33}` components over finite Toy I
molecules with $|M_U|=n$ and all permitted tie-breaks. Three explicit
families attain equality for every $n\ge4$. The paper also records the
scope-limited consequence

$$
v_d(M)\ge\left\lceil\frac{\rho(M)}{3}\right\rceil-10d
$$

within Toy I.

## Proof method

The all-size result is a combinatorial amortised-flow proof, not an inference
from finite enumeration.

1. Exact event accounting gives $T=n-1-U$, where $T$ is the number of `{33}`
   components and $U$ is the number of upper-layer atoms whose cross-layer
   partner was removed earlier.
2. For each lower-layer cleanup, a local parent-slot count bounds its orphaned
   cross-bonds. Bottom fixed ends created by earlier cuts supply injectively
   matched tokens, yielding $U\le2T$.
3. Combining the two identities gives the lower bound. Three explicit
   molecule families give equality, proving sharpness.

The result does **not** determine the intrinsic optimum $W$, classify its
complexity, cover arbitrary legal cut sequences, or extend to Toy I plus/II/III
or the full molecule problem.

## Contents

- `MANUSCRIPT.pdf` — seven-page preprint; `MANUSCRIPT.tex` is its canonical source;
- `MANUSCRIPT.md` — line-oriented reading copy;
- `src/molecule_cut/` — model, cutting algorithm, enumerators, and exact DP;
- `tests/` — unit, family, exhaustive-small-slice, and regression tests;
- `scripts/` and `data/` — deterministic replay scripts and stored tables;
- `REPRODUCIBILITY.md` — verification commands and finite-audit scope;
- `CITATION.cff` and `LICENSE` — citation metadata and code/data licensing.

`e-jc.sty` is the unchanged official E-JC style file used to compile the paper.
The manuscript itself remains under the author's copyright; see `LICENSE` for
the precise scope of the MIT code/data license.

## Reproduce

Requires Python 3.12 and [uv](https://docs.astral.sh/uv/).

```sh
uv sync --frozen --group dev
uv run --frozen pytest
uv run --frozen ruff check .
```

The tests are an executable regression layer.  The all-size result is proved in
the manuscript; it is not inferred from finite enumeration.

## Source and citation

The source model is Deng--Hani--Ma, arXiv:2408.07818v3, Sections 11.1--11.2:
<https://arxiv.org/abs/2408.07818>.

If this code or preprint is useful, please cite the release metadata in
`CITATION.cff`.  The target-journal author guidelines for the Electronic Journal
of Combinatorics are at
<https://www.combinatorics.org/ojs/index.php/eljc/about/submissions#authorGuidelines>.
