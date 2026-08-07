# molecule-cut

Reproducible code and manuscript for the sharp Toy I benchmark in the
Deng--Hani--Ma two-layer molecule cutting problem.

The paper proves, for every finite Toy I molecule and every permitted execution
of the source Definition 11.4 algorithm,

\[
\#\{33\}\ge \left\lceil\frac{|M_U|-1}{3}\right\rceil.
\]

It also gives explicit families attaining the bound and records the resulting
scope-limited intrinsic \(v_d\) consequence. The manuscript does not claim a
global theorem for the full molecule, Toy I plus/II/III, arbitrary legal cuts,
or complexity hardness.

## Contents

- `MANUSCRIPT.pdf` — compiled paper;
- `MANUSCRIPT.tex` and `MANUSCRIPT.md` — paper sources;
- `SPEC.md` — public implementation conventions;
- `e-jc.sty` — unchanged official E-JC style file used for the submission PDF;
- `EJC-SUBMISSION-STATEMENT.md` — originality and submission checklist;
- `src/molecule_cut/` — model, cutting algorithm, enumerators, and exact DP;
- `tests/` — unit, exhaustive-small-slice, family, and regression tests;
- `scripts/` and `data/` — bounded reproduction scripts and stored tables;
- `REPRODUCIBILITY.md` — commands and finite-audit scope.

`MANUSCRIPT.tex` is the canonical source for the compiled PDF; `MANUSCRIPT.md`
is the line-oriented reading copy.  Their displayed theorem numbers can differ
because the E-JC style assigns a shared counter to lemma/theorem environments.

The source model is Deng--Hani--Ma, arXiv:2408.07818v3, Sections 11.1--11.2:
<https://arxiv.org/abs/2408.07818>.

The target-journal author guidelines are available at
<https://www.combinatorics.org/ojs/index.php/eljc/about/submissions#authorGuidelines>.

## Quick start

Requires Python 3.12 and `uv`.

```sh
uv sync --frozen
uv run --frozen pytest
uv run --frozen ruff check .
```

The tests are the executable regression layer; the all-size theorem is proved
in the manuscript and is not inferred from finite enumeration.

The current release candidate includes the upper-layer degree lemma and its
finite guard.  The full regression suite reports 653 passing tests; the
focused `test_p_exact_route.py` audit reports 4 passing tests.
