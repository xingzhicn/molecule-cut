# Public code specification

This file gives the conventions used by the public implementation. The source
model is Deng--Hani--Ma, arXiv:2408.07818v3, Sections 11.1--11.2.

## 1. Molecules

A molecule is a directed graph whose vertices are atoms and whose edges are
bonds. Each atom has four slots. Degree is the number of incident bonds plus
free ends; fixed ends do not count toward degree. Each atom has at most two
parents and at most two children. Cutting a connected set as free turns each
boundary bond into one persistent fixed end at the surviving endpoint.

In Toy I, the upper and lower layers are trees with no initial fixed ends, no
lower atom is a parent of an upper atom, and every upper atom has exactly one
cross-layer bond directed to the lower layer.

## 2. Elementary pieces and score

The implementation recognises the elementary singleton and pair types used by
the source algorithm: `{4}`, `{3}`, `{33}`, and the zero-score two-fixed-end
singleton. A `{343}` pattern is a cleanup macro that is split into one `{3}`
and one `{33}`. For a complete cut `C`,

```text
score_d(C) = #{33}(C) - 10 d #{4}(C).
```

The circuit rank of a connected graph is `rho = |E| - |V| + 1`.

## 3. Algorithms and quantities

The source six-step Toy I algorithm is implemented by
`molecule_cut.algorithm.run_toy1_algorithm`, including its permitted
tie-breaks. The exact intrinsic dynamic program is exposed by
`molecule_cut.exhaustive.v_d` and `w`.

The manuscript distinguishes the benchmark minimum over source-algorithm
executions from the intrinsic maximum over all legal complete cuts. Finite
enumeration is used only as regression evidence; the all-size `1/3` theorem is
the written amortised-flow proof in `MANUSCRIPT.pdf`.
