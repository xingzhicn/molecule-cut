# The sharp benchmark constant for the Deng--Hani--Ma two-layer cutting algorithm

**Status.** Preprint release, 2026-08-24. The benchmark theorem below
is proved for the finite Toy I class and every legal tie-break of Definition
11.4. A separate section then records the scope-limited consequence for the
intrinsic optimum `v_d`. Finite computation is regression evidence only.

## Abstract

The simplified two-layer molecule argument of Deng--Hani--Ma gives
`#{33} >= (|M_U|-1)/5` for its Toy I cutting algorithm. We isolate the exact
source of the loss and prove the stronger, sharp benchmark bound

```text
#{33} >= ceil((|M_U|-1)/3).
```

The proof is an amortised flow argument. A cleanup cut can orphan cross-bonds,
and an operation (b) can orphan one additional cross-bond that is invisible to
the older bookkeeping. For events occurring while an upper-layer atom remains,
local slot capacity is paid for by bottom-fixed-end tokens created by earlier
cuts. An injective token map cancels the source terms and gives
`#CUT_U <= 2#{33}`. Together with the exact identity
`#{33}=|M_U|-1-#CUT_U`, this proves the lower bound. Three explicit families give
the matching upper bound, hence

```text
a_n = ceil((n-1)/3),  n >= 4.
```

This is a combinatorial benchmark result. It neither improves the Boltzmann
theorem itself nor asserts a global no-go. As a direct consequence of the
benchmark witness, every finite Toy I molecule also satisfies the scoped
intrinsic lower bound

```text
v_d(M) >= ceil(rho(M)/3) - 10d,
```

with the asymptotic interpretation stated in Section 4.

## 1. Source contract and scope

We use the two-layer model and Toy I in Section 11 of the public arXiv version
of Deng--Hani--Ma [DHM, Definitions 11.1--11.4]. Each atom has four slots,
partitioned into bonds, free ends, and fixed ends; fixed ends do not count
toward degree. There are at most two parent slots and at most two child slots.
A cross-layer bond is directed from `M_U` to `M_D` and occupies a parent slot of
its lower endpoint. Cutting a connected set as free turns each boundary bond
into one persistent fixed end at the surviving endpoint (bottom at a surviving
parent and top at a surviving child); fixed ends are neither rewired nor
duplicated. Both layers are trees with no initial fixed ends; no atom of `M_D`
is a parent of an atom of `M_U`; and every atom of `M_U` has exactly one bond to
`M_D`. The degree convention is

```text
degree = number of bonds + number of free ends;
fixed ends do not count toward degree.
```

The residual `M_D` is called proper as in [DHM, Definition 11.3]. We run the
six-step algorithm of [DHM, Definition 11.4], allowing every permitted choice.
Proposition 11.7 guarantees that a complete run cuts only elementary
components and terminates.

Call a complete execution of the six-step algorithm, with any permitted
tie-break and an empty final residual molecule, a **Definition-11.4 execution**.
Its output is elementary after the prescribed `{343}` cleanup macro is resolved
into one `{3}` and one `{33}`; the macro contributes one to `#{33}` and no
`#{4}`. Separately, a **legal complete cut** means any sequence of elementary
cuts that ends with an empty
residual molecule. The intrinsic quantity `v_d` maximizes the resulting score
over these general legal complete cuts, whereas the benchmark theorem concerns
only Definition-11.4 executions.

For a run, let `T=#{33}` and `n=|M_U|`. Define the benchmark

```text
a_n = min over finite Toy I molecules with |M_U|=n,
      min over legal Definition-11.4 tie-breaks, of T.
```

The source result is `#{4}=1` and `T >= (n-1)/5`. Our result replaces only the
second estimate for this benchmark. The quantities

```text
v_d(M) = max over all legal complete cuts (#{33}-10d·#{4}),
W(M)   = max over all legal complete cuts #{33}
```

are separate project objects. The benchmark theorem addresses `a_n`; the
scope-limited transfer to `v_d` is proved in Section 4. It does not determine
`W`, other toy models, or arbitrary legal cutting sequences.

### 1.1 Exact public source location

The source definitions and the original `1/5` estimate are in the public
arXiv version of Deng--Hani--Ma, Sections 11.1--11.2, especially Definitions
11.1--11.4 and Proposition 11.7:
<https://arxiv.org/abs/2408.07818>. The present paper restates every graph,
slot, cut, and quantifier needed for its new counting argument.

### 1.2 Relation to prior work

Deng--Hani--Ma supply the original `1/5` benchmark estimate and its role in
the hard-sphere derivation; the survey of Bodineau--Gallagher--Saint-Raymond--
Simonella reproduces the surrounding reduction. Bruned--Clarisse study a
Kruskal-style reduction for the wave-molecule problem, not the quantitative
hard-sphere count proved here. Our claim is therefore a quantitative Toy I
benchmark result, not an algorithm-classification claim or a new kinetic
limit theorem. The literature search was bounded; we make no exhaustive
priority claim.

## 2. Exact accounting

Call an operation `(a)` an `A` event, operation (b) a `B` event, a cleanup of
an adjacent pair a `C` event (`CUT33`), a cleanup of a `{343}` triple a `D` event
(`CUT343`), and a one-atom upper-layer removal after its cross-bond partner has
already gone a `U` event (`CUT_U`). Let `E_<` be the subset of `B`, `C`, and `D`
events whose pre-cut state still contains at least one atom of `M_U`. This
cutoff is part of the statement: Step (6) cleanups after the upper layer is
empty cannot orphan a future upper-layer partner.

The first identity in (2.0) below asserts that every `B` event contributes
exactly one `{33}`. That requires the selected upper atom to have degree three,
which the following lemma supplies. The source asserts the conclusion — in the
proof of [DHM, Proposition 11.7], each cut atom of `M_U` "either forms a
component by itself or belongs to a `{33}` molecule" — but its step (2) also
carries the defensive clause "cut `m` as free from `{n,m}` if `n` has deg 4".
The two are consistent only because that clause is vacuous on Toy I, which we
now prove.

### Lemma 0 (upper-layer degree at operation (b))

In a Definition-11.4 execution on a finite Toy I molecule, let `n in M_U` be
the atom selected at step (1) and suppose its cross-bond partner `m in M_D` has
not yet been cut. Then `n` has degree three, unless `n` is the first atom
selected by the whole run, in which case `m` has degree four. Consequently
every `B` event cuts a genuine `{33}`, and the deg-4 clause of step (2) is
never invoked.

**Proof.** In Toy I every atom of `M_U` carries exactly one cross-bond, which
occupies one of its at most two child slots; hence every atom of `M_U` has at
most one child in `M_U`. Summing the number of `M_U`-children over `M_U`
counts the `|M_U|-1` intralayer bonds, so with every term at most one, exactly
one atom `n_0 in M_U` has no `M_U`-child.

No top fixed end ever appears in `M_U`: step (1) selects only atoms whose
`M_U`-children have all been cut, so cutting an upper atom creates no top fixed
end above it, and cuts inside `M_D` create none either. Hence the fixed ends at
the selected `n` are exactly the bottom fixed ends created by cutting children
of `n`. By hypothesis `m` survives, so the cross-bond contributes none, and
the remaining children of `n` lie in `M_U` and have all been cut. Therefore

```text
deg(n) = 4 - #{M_U-children of n} = 3 if n != n_0, and 4 if n = n_0.
```

At the start nothing has been cut, so `n_0` is the unique atom satisfying the
selection rule of step (1); it is the first atom the run selects. At that
moment `M_D` is still a tree with no fixed ends, so every atom of `M_D` has
degree four. Hence `m` has degree four and step (2) takes its "otherwise"
branch, which is an `A` event. Every selection with a degree-three partner
therefore has `n != n_0` and `deg(n)=3`. ∎

For a complete run, the source counting identities used below are

```text
T = B + C + D,       A - (C + D) = 1,       A + B + U = n,       (2.0)
```

where `T=#{33}`. The first identity counts the `{33}` produced by `(b)`,
`CUT33`, and the `{343}` split; the second is the change of the source
potential `lambda`; the third partitions the upper-layer atoms. The first two
are the source's `{33}` and `lambda` identities, while the third is the elementary
upper-layer partition; all are stated here so the subsequent elimination is
self-contained.

For `e in E_<`, let `S_e` be its `M_D` cut set. Define:

* `q_e` as the number of live cross-bonds removed by `e`, except that for a `B`
  event the selected `U--D` bond is omitted;
* `s_e` as the number of live `M_D` boundary arcs `p -> v` with `p` outside
  `S_e` and `v` in `S_e`;
* `r_e` as the number of degree-three atoms in `S_e` carrying a bottom fixed end.

Thus `q_e` is the number of future `CUT_U` atoms orphaned by the event, while
each `s_e` creates one bottom-fixed-end token at the surviving `M_D` parent.

### Lemma 1 (exact orphan partition)

For every complete run,

```text
#CUT_U = sum_{e in E_<} q_e.                              (2.1)
```

**Proof.** Every `CUT_U` atom has lost its unique initial cross-bond partner
earlier. A `C` or `D` event removes exactly the live cross-bonds counted by
`q_e`. A `B` event removes its selected upper atom together with the selected
bond, and its other live cross-bonds are exactly those counted by `q_e`.
Conversely, each counted bond has one still-present upper endpoint and produces
one later `CUT_U`. Events after the upper layer is empty have no live cross-bond
to a future `CUT_U`, proving the partition. ∎

### Lemma 2 (local flow capacity)

For every `e in E_<`,

```text
q_e + s_e <= 2 + r_e.                                      (2.2)
```

**Proof.** For a degree-three `M_D` atom `v`, write `c_v` for its live
cross-bonds, `p_v` for its live `M_D` parents, and `t_v` for `1` or `0`
according as its unique fixed end is bottom or top. The atom has two parent
slots; a top fixed end occupies one of them, while a bottom fixed end does not.
Thus its two parent slots imply

```text
c_v + p_v + (1-t_v) <= 2.                                  (2.3)
```

For `B`, the selected cross-bond is one of the `c_v`, so
`q_e+s_e <= t_v`, which is stronger than (2.2). For `C`, the one internal
`M_D` bond contributes one parent slot across the two atoms; summing (2.3)
gives `q_e+s_e <= 1+r_e`. For `D`, the two outer atoms are degree three and the
centre is degree four with no fixed end, so its two parent slots give
`c_v+p_v<=2`. The two internal bonds contribute two parent slots. Summing the
three slot inequalities gives
`q_e+s_e <= 2+r_e`. These cases cover all orientations and allow double-parent
atoms in either layer. ∎

### Lemma 3 (bottom-token injection)

```text
sum_{e in E_<} r_e <= sum_{e in E_<} s_e.                    (2.4)
```

**Proof.** Initially `M_D` has no fixed ends. Before an event in `E_<`, the
only operations that have removed `M_D` atoms are earlier `B`, `C`, or `D`
events. An operation (a) or `CUT_U` removes only an upper atom and creates a
top fixed end on the residual lower side; `CUT3` is Step (6) and cannot occur
while an upper atom remains. Hence every bottom fixed end counted by an `r_e`
has a unique earlier source boundary arc counted by some `s_f`. The cut removes
that bond and never duplicates or rewires its fixed end, so distinct consumed
tokens have distinct source arcs. ∎

### Theorem 4 (sharp Toy I benchmark bound)

For every finite Toy I molecule and every legal tie-break of Definition 11.4,

```text
#{33} >= ceil((|M_U|-1)/3).                                  (2.5)
```

Consequently, `a_n >= ceil((n-1)/3)`.

**Proof.** Sum (2.2) over `E_<` and use (2.4):

```text
sum q_e + sum s_e <= 2|E_<| + sum r_e
                 <= 2|E_<| + sum s_e.
```

After cancellation, (2.1) gives `#CUT_U <= 2|E_<|`. Every `B`, `C`, or `D`
event in `E_<` produces one `{33}` component, so `|E_<| <= #{33}` and

```text
#CUT_U <= 2#{33}.                                            (2.6)
```

Let `A,B,C,D,U` count the events defined above. Substituting (2.0) and
eliminating `A,B,C,D` gives the exact project identity

```text
#{33} = |M_U| - 1 - #CUT_U.                                  (2.7)
```

Combining (2.6) and (2.7) yields `|M_U|-1 <= 3#{33}`. Since `#{33}` is an
integer, (2.5) follows. ∎

### 2.1 The exact contract and the old route

Writing `b_k` for the cross-bonds broken by cleanup `k`, and `q` for the sum
of the additional `B`-event orphans, (2.1) is

```text
#CUT_U = sum_k b_k + q.
```

Thus the exact sufficient-and-necessary run-level contract for (2.5) is

```text
P*: sum_k b_k + q <= 2(B + C + D).                           (2.8)
```

The older proposed condition `P: sum_k b_k <= 2(C+D)` is sufficient, but it is
not equivalent to (2.8), and it is not used or claimed here. The `q` correction
is real: a bounded witness in `tests/test_p_exact_route.py` has a `B` event with
`q=1`.

## 3. Sharpness: matching explicit families

The following families are finite Toy I molecules, defined directly here.

**Transpose family.** For `k>=1`, take `n=3k+1`, let
`M_U=(u_0,...,u_{3k})` and `M_D=(d_0,...,d_{3k-1})` be directed paths with
internal bonds `u_i->u_{i+1}` and `d_j->d_{j+1}`, and write
`i=qk+r` (`q=0,1,2`, `0<=r<k`) for `0<=i<3k`. Join `u_i` to `d_{3r+q}`,
and add `u_{3k}->d_0`. The lowest upper atoms are removed in the order
`u_{3k},u_{3k-1},...,u_{2k}`, giving `k+1` operations `(a)`. Their lower
partners are `d_0,d_{3r+2}`. After these cuts the degree-three set is `{d_0}`
together with `{d_{3r+2}:0<=r<k}`; the points are separated along the path
until the final cut, when `d_0,d_2` and their degree-four neighbour `d_1` form
the first `{343}`. The blocks `{d_{3r},d_{3r+1},d_{3r+2}}` can then be removed
successively in the indicated legal tie-break: after block `r` is removed, `d_{3r+3}` becomes degree three,
`d_{3r+5}` was already degree three, and `d_{3r+4}` remains degree four because
its cross-bond comes from the still present `u_{k+r+1}`. For `0<=r<k-1`,
these facts give the next block; the final block `r=k-1` is terminal and needs
no successor indices. In block `r`, the three designated cross-bonds come from
`u_r,u_{k+r},u_{2k+r}`. For `r>=1`, only the last endpoint has already been
removed, so two live cross-bonds are broken. For `r=0`, the extra bond
`u_{3k}->d_0` is also present, but both it and the designated `u_{2k}` endpoint
have already been removed, so again exactly two live cross-bonds are broken.
This gives the block cascade by induction. After the `k` cleanups, the remaining
`2k` upper atoms are `CUT_U`, and the profile is
`(A,B,C,D,CUT_U)=(k+1,0,0,k,2k)`, giving `#{33}=k` and `#{4}=1`.

**Padded family.** Append one atom to the top of `M_U`, one to the end of
`M_D`, and join the new atoms. The old lower endpoint trades a free end for a
child and keeps degree four. The new upper endpoint is cut last, so the
transpose phase and its `k` cleanups are unchanged. The final lower atom then
has degree three and closes with one operation `(b)`. For `n=3k+2`, the
profile is `(A,B,C,D,CUT_U)=(k+1,1,0,k,2k)` and `#{33}=k+1`.

**Trimmed family.** For `n=3k` with `k>=2`, take paths
`M_U=(u_0,...,u_{3k-1})` and `M_D=(d_0,...,d_{3k-2})`. Use the transpose
cross-bonds except redirect the bond that would end at `d_{3k-1}` to `d_0`;
omit `u_{3k}`. The first `k` upper cuts and next `k-1` block cleanups are as
above. The cascade stops at the tail pair `{d_{3k-3},d_{3k-2}}`; one further
operation `(a)` lowers the latter to degree three, and the pair is then a
`{33}` cleanup. This final cleanup breaks one live cross-bond, while each prior
block breaks two. The remaining `2k-1` upper atoms are `CUT_U`, so the profile
is `(A,B,C,D,CUT_U)=(k+1,0,1,k-1,2k-1)` and `#{33}=k`.

Therefore `a_n <= ceil((n-1)/3)` for every `n>=4`. Together with Theorem 4:

### Corollary 5

```text
a_n = ceil((n-1)/3),  n >= 4.                                (3.1)
```

This is sharp for the Definition-11.4 benchmark, including its worst molecule
and worst tie-break quantifiers. It is not a statement about the intrinsic
optimum over all cutting sequences.

## 4. Intrinsic consequence on Toy I

Theorem 4 is a statement about a prescribed algorithm. The following
proposition uses that algorithm only as one legal witness in the intrinsic
maximum; it does not assert that every legal cutting sequence has the same
count.

### Proposition 6 (scoped intrinsic lower bound)

Fix `d > 0`. Let `M` be a finite Toy I molecule with both layers nonempty and
write `n=|M_U|`. Then

```text
v_d(M) >= ceil((n-1)/3) - 10d
        = ceil(rho(M)/3) - 10d.                           (4.1)
```

**Proof.** Put `m=|M_D|`. Each layer is a tree, so the two intralayer edge
counts are `n-1` and `m-1`. Toy I gives exactly one cross-bond for each upper
atom, hence exactly `n` cross-bonds. The full graph is connected. Here `rho(M)`
denotes its circuit rank, `rho(M)=|E|-|V|+1`, and therefore

```text
rho(M) = |E|-|V|+1
        = ((n-1)+(m-1)+n)-(n+m)+1
        = n-1.                                             (4.2)
```

Proposition 11.7 of [DHM] supplies a legal complete Definition-11.4 run with
`#{4}=1`. Theorem 4 gives `#{33} >= ceil((n-1)/3)` for that run.
Since `v_d` maximizes `#{33}-10d·#{4}` over all legal complete cuts, it is at
least the score of this one run, which proves (4.1). ∎

### Corollary 7 (asymptotic Toy I class bound)

For `n>=4`, define the finite-size intrinsic constant over the complete Toy I
class by an infimum (the class has no a priori bound on `|M_D|`):

```text
c_d^ToyI(n) = inf over finite Toy I M with |M_U|=n of v_d(M)/rho(M).
```

Then

```text
c_d^ToyI(n) >=
  (ceil((n-1)/3)-10d)/(n-1).                              (4.3)
```

Consequently, for fixed `d`,

```text
liminf_{n->infinity} c_d^ToyI(n) >= 1/3.                  (4.4)
```

This is an asymptotic statement within Toy I. It is not the global infimum
`c_d^*`, whose definition includes small positive circuit ranks and whose
positive lower bound is not claimed here.

## 5. Adversarial audit and reproducibility

An independent adversarial read checks the exact places where the old route
could silently overclaim:

1. the `E_<` cutoff excludes Step (6) token creation;
2. the extra `B`-event orphan `q` is retained;
3. the local slot count uses parent slots and works for all three orientations
   of `{343}` and for double-parent configurations;
4. bottom tokens are injected only from earlier lower-layer removals;
5. finite enumeration is used only as regression, never as the infinite proof.

The bounded full-labelled audit covers 11,590 molecules and 25,243 complete
tie-break records with `|M_U|+|M_D|<=6`; it includes double-parent cases in both
layers and reports no failed records or `P*` violations. These numbers are
finite evidence only. The public code and reproduction repository is
<https://github.com/xingzhicn/molecule-cut>. The commands are listed below.

```text
uv run --frozen pytest
uv run --frozen ruff check .
```

## 6. What is and is not concluded

Concluded:

* the source's Toy I/Definition-11.4 benchmark has the exact constant `1/3`;
* the proof is an all-size combinatorial argument, independent of the finite
  audit and of the single-parent/chain restrictions used by older experiments;
* the old `P` route is unnecessary, while `P*` is the exact accounting contract;
* the intrinsic consequence (4.1) and the asymptotic class bound (4.4) hold
  under the explicitly stated Toy I scope.

Not concluded:

* no global ``no-localization'', ``no-beam'', or ``B-v0.1 impossible'' claim;
* no global positive lower bound for `c_d^*`;
* no exact formula, complexity classification, or optimal-cut theorem for `W`,
  Toy I plus/II/III, or the full molecule;
* no claim that the score bound holds for every arbitrary legal cutting
  sequence; Proposition 6 uses one prescribed legal witness;
* no complexity classification;
* no claim that the finite audit proves an infinite statement.

## References

* **[DHM]** Y. Deng, Z. Hani, X. Ma, *Long time derivation of the Boltzmann
  equation from hard sphere dynamics*, arXiv:2408.07818v3 (2025), Sections
  11.1--11.2. <https://arxiv.org/abs/2408.07818>.
* **[BGSS]** T. Bodineau, I. Gallagher, L. Saint-Raymond, S. Simonella,
  *Derivation of the Boltzmann equation from hard-sphere dynamics (after Deng,
  Hani, Ma)*, Séminaire Bourbaki 1247, arXiv:2602.04407.
* **[BC]** Y. Bruned, V. Clarisse, *Kruskal-style algorithm for cubic
  Schrödinger equation molecule reduction*, arXiv:2603.23298.
