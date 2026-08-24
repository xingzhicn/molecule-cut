"""Finite full-domain guards for the exact `1/3` reduction.

This is deliberately a bounded computation, not a proof of P or P*.  Unlike
the older random checks, it ranges over both double-parent orientations and
every tie-break branch in the labelled total-size-six slice.
"""

from molecule_cut.algorithm import OpKind
from molecule_cut.builders import build
from molecule_cut.enumerate import enumerate_toy1_full_labeled
from molecule_cut.exhaustive import enumerate_tiebreak_records
from molecule_cut.fast_vd import fast_w
from molecule_cut.molecule import Direction, Layer

CLEANUP_OPS = frozenset({OpKind.CUT33, OpKind.CUT343})
FLOW_OPS = CLEANUP_OPS | {OpKind.B}


def _assert_toy1_contract(mol, n_up: int) -> tuple[bool, bool]:
    """Check the structural domain conditions used by the finite audit."""
    mol.check_invariants()
    up = set(mol.layer_atoms(Layer.UP))
    down = set(mol.layer_atoms(Layer.DOWN))
    assert len(mol.cross_bonds()) == n_up
    assert mol.circuit_rank() == n_up - 1
    assert len(mol.components(within=up)) == 1
    assert len(mol.components(within=down)) == 1
    assert all(atom.n_fixed() == 0 for atom in mol.atoms.values())
    assert all(sum(child in down for child in mol.children(u)) == 1 for u in up)

    up_double_parent = any(sum(parent in up for parent in mol.parents(a)) == 2 for a in up)
    down_double_parent = any(
        sum(parent in down for parent in mol.parents(a)) == 2 for a in down
    )
    return up_double_parent, down_double_parent


def _reduction_statistics(rec) -> tuple[int, int, int, int, int]:
    """Return (cleanup_bonds, B_extra_bonds, cleanups, B_ops, CUT_U)."""
    cleanup_bonds = cleanups = b_extra_bonds = b_ops = cut_u = 0
    for step in rec.steps:
        if step.op in CLEANUP_OPS:
            cleanups += 1
            cleanup_bonds += step.cross_bonds_broken
        elif step.op is OpKind.B:
            # The chosen U--D cross bond lies inside this B cut; any other
            # broken live cross bond leaves a future U atom partnerless.
            assert step.cross_bonds_broken >= 1
            b_ops += 1
            b_extra_bonds += step.cross_bonds_broken - 1
        elif step.op is OpKind.CUT_U_ONLY:
            cut_u += 1
    return cleanup_bonds, b_extra_bonds, cleanups, b_ops, cut_u


def _bottom_fixed_support(state, subset: frozenset[int]) -> int:
    """Count degree-three MD atoms in a cleanup with a bottom fixed end."""
    return sum(
        state.atoms[atom].layer is Layer.DOWN
        and state.degree(atom) == 3
        and Direction.BOTTOM in state.atoms[atom].fixed_directions()
        for atom in subset
    )


def _flow_event_statistics(state, step) -> tuple[int, int, int]:
    """Return (orphan_crosses, new_bottom_tokens, consumed_bottom_tokens)."""
    subset = set(step.subset)
    down = set(state.layer_atoms(Layer.DOWN))
    if step.op is OpKind.B:
        orphan_crosses = step.cross_bonds_broken - 1
    else:
        orphan_crosses = step.cross_bonds_broken
    new_bottom_tokens = sum(
        bond.parent not in subset
        and bond.child in subset
        and bond.parent in down
        and bond.child in down
        for bond in state.bonds
    )
    consumed_bottom_tokens = _bottom_fixed_support(state, step.subset)
    return orphan_crosses, new_bottom_tokens, consumed_bottom_tokens


def test_exact_reduction_on_all_labelled_toy1_runs_of_total_size_at_most_six():
    """[Exact computation] Full labelled/all-tie audit on the finite slice T_6."""
    molecules = records = 0
    max_branches = 0
    up_double_parent = down_double_parent = False

    for n_up in range(1, 6):
        for n_down in range(1, 6):
            if n_up + n_down > 6:
                continue
            for mol in enumerate_toy1_full_labeled(n_up, n_down):
                molecules += 1
                has_up_double, has_down_double = _assert_toy1_contract(mol, n_up)
                up_double_parent |= has_up_double
                down_double_parent |= has_down_double

                branches = enumerate_tiebreak_records(mol)
                max_branches = max(max_branches, len(branches))
                for _, rec in branches:
                    records += 1
                    assert rec.failed is None, rec.failed
                    cleanup_bonds, b_extra_bonds, cleanups, b_ops, cut_u = (
                        _reduction_statistics(rec)
                    )

                    # Lemma 1 and the two proved counting identities.
                    assert cut_u == cleanup_bonds + b_extra_bonds
                    assert rec.n33 == cleanups + b_ops

                    # P is only sufficient; P* is the exact target.  Both
                    # hold on this finite, fully enumerated slice.
                    assert cleanup_bonds <= 2 * cleanups
                    assert cleanup_bonds + b_extra_bonds <= 2 * (cleanups + b_ops)

                    # Full finite replay of the flow proof.  Only events whose
                    # pre-cut state still contains U atoms enter its token
                    # injection; later cleanups have no live cross bonds.
                    state = mol.copy()
                    flow_orphans = flow_sources = flow_consumed = flow_events = 0
                    for step in rec.steps:
                        if step.op is OpKind.CUT3:
                            assert not state.layer_atoms(Layer.UP)
                        if step.op in CLEANUP_OPS:
                            support = _bottom_fixed_support(state, step.subset)
                            assert max(step.cross_bonds_broken - 2, 0) <= support
                        if step.op in FLOW_OPS:
                            q, sources, consumed = _flow_event_statistics(state, step)
                            if state.layer_atoms(Layer.UP):
                                assert q + sources <= 2 + consumed
                                flow_orphans += q
                                flow_sources += sources
                                flow_consumed += consumed
                                flow_events += 1
                                assert flow_consumed <= flow_sources
                            else:
                                assert q == 0
                        state = state.cut_as_free(set(step.subset))

                    assert flow_orphans == cut_u
                    assert flow_orphans <= 2 * flow_events
                    assert flow_events <= rec.n33

    # Fixed counts prevent either the domain or tie-break traversal from being
    # silently narrowed.  They are not extrapolated beyond this finite slice.
    assert molecules == 11_590
    assert records == 25_243
    assert max_branches == 9
    assert up_double_parent
    assert down_double_parent


UPPER_OPS = frozenset({OpKind.A, OpKind.B, OpKind.CUT_U_ONLY})


def test_lemma_0_upper_atom_has_degree_three_at_every_b_event():
    """[Exact computation] Finite guard for Lemma 0 of MANUSCRIPT section 2.

    Definition 11.4 step (2) carries a "cut m as free from {n,m} if n has deg 4"
    clause.  Lemma 0 proves the clause is vacuous on Toy I: the selected upper
    atom has degree four only for the very first selection of a run, and that
    selection always meets a degree-four partner, so it is an A event.  If the
    clause were reachable, a B cut would not be a {33} and the identity
    T = B + C + D behind (2.0) would fail.
    """
    b_events = deg_four_selections = 0

    for n_up in range(1, 6):
        for n_down in range(1, 6):
            if n_up + n_down > 6:
                continue
            for mol in enumerate_toy1_full_labeled(n_up, n_down):
                for _, rec in enumerate_tiebreak_records(mol):
                    assert rec.failed is None, rec.failed
                    state = mol.copy()
                    for index, step in enumerate(rec.steps):
                        if step.op in UPPER_OPS:
                            (upper,) = [
                                atom
                                for atom in step.subset
                                if state.atoms[atom].layer is Layer.UP
                            ]
                            degree = state.degree(upper)
                            if step.op is OpKind.B:
                                b_events += 1
                                assert degree == 3, (
                                    f"B event with deg(n)={degree}: the cut is not a "
                                    f"{{33}} and Lemma 0 fails"
                                )
                            if degree == 4:
                                deg_four_selections += 1
                                assert index == 0, (
                                    f"deg-4 upper atom selected at step {index}, "
                                    f"not only at the first selection"
                                )
                                assert step.op is OpKind.A, step.op
                        state = state.cut_as_free(set(step.subset))

    # Every run selects exactly one degree-four upper atom: its first.
    assert deg_four_selections == 25_243
    assert b_events == 1_530


def test_full_toy1_tie_breaking_can_change_good_component_and_cleanup_counts():
    """A full-domain witness rules out global tie-invariance claims.

    The intrinsic optimum is still well-defined and equals 3 here; only the
    Def. 1.4 greedy execution depends on its cleanup tie-break.
    """
    mol = build(
        4,
        3,
        [(1, 0), (2, 0), (0, 3)],
        [(1, 0), (0, 2)],
        [(0, 1), (1, 0), (2, 2), (3, 1)],
    )
    outcomes = set()
    for _, rec in enumerate_tiebreak_records(mol):
        assert rec.failed is None, rec.failed
        cleanup_bonds, _, cleanups, _, _ = _reduction_statistics(rec)
        outcomes.add((rec.n33, cleanups, cleanup_bonds))

    assert outcomes == {(2, 1, 1), (3, 0, 0)}
    assert fast_w(mol) == 3


def _full_domain_tie_break_summary(n_up: int, n_down: int):
    """Return exact labelled counts for a small full-domain tie-break audit."""
    molecules = multi = varied = 0
    outcome_counts = {}
    for mol in enumerate_toy1_full_labeled(n_up, n_down):
        molecules += 1
        records = enumerate_tiebreak_records(mol)
        assert all(record.failed is None for _, record in records)
        values = tuple(sorted({record.n33 for _, record in records}))
        outcome_counts[values] = outcome_counts.get(values, 0) + 1
        multi += len(records) > 1
        varied += len(values) > 1
    return molecules, multi, varied, outcome_counts


def test_full_toy1_n3_exactly_refutes_f12_invariance():
    """All labelled (3,3) Toy I instances include many non-invariant branches."""
    assert _full_domain_tie_break_summary(3, 3) == (
        1_134,
        378,
        252,
        {(1,): 468, (2,): 414, (1, 2): 252},
    )


def test_full_toy1_n4_n3_exactly_refutes_f12_invariance():
    """The same failure persists on the complete labelled (4,3) slice."""
    assert _full_domain_tie_break_summary(4, 3) == (
        7_560,
        4_536,
        2_112,
        {
            (1,): 1_800,
            (2,): 3_288,
            (3,): 360,
            (1, 2): 1_320,
            (1, 3): 48,
            (2, 3): 672,
            (1, 2, 3): 72,
        },
    )


def test_b_operation_correction_is_a_real_part_of_the_exact_reduction():
    """A B step can orphan one additional U atom, so q cannot be dropped."""
    mol = build(
        4,
        3,
        [(1, 0), (2, 0), (0, 3)],
        [(1, 0), (0, 2)],
        [(0, 0), (1, 1), (2, 1), (3, 2)],
    )
    outcomes = set()
    for _, rec in enumerate_tiebreak_records(mol):
        assert rec.failed is None, rec.failed
        cleanup_bonds, b_extra_bonds, cleanups, b_ops, cut_u = _reduction_statistics(rec)
        outcomes.add((cleanup_bonds, b_extra_bonds, cleanups, b_ops, cut_u, rec.n33))

    assert (0, 1, 1, 1, 1, 2) in outcomes
