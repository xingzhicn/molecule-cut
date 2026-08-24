"""Finite replays for the G1 block-orientation encoding lemma."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations, product

import pytest

from molecule_cut.algorithm import run_toy1_algorithm
from molecule_cut.builders import build
from molecule_cut.elementary import ElemType, classify
from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.exhaustive import _legal_moves
from molecule_cut.families import padded_family, transpose_family, trimmed_family
from molecule_cut.molecule import Atom, Bond, Direction, EndKind, Molecule


def _complete_sequences(mol: Molecule) -> list[list[tuple[frozenset[int], ElemType]]]:
    if not mol.atoms:
        return [[]]
    sequences: list[list[tuple[frozenset[int], ElemType]]] = []
    for block in _legal_moves(mol):
        kind = classify(mol, set(block))
        for suffix in _complete_sequences(mol.cut_as_free(set(block))):
            sequences.append([(block, kind), *suffix])
    return sequences


def _quotient_order(sequence: list[tuple[frozenset[int], ElemType]], mol: Molecule) -> list[int]:
    owner = {atom: index for index, (block, _) in enumerate(sequence) for atom in block}
    successors: dict[int, set[int]] = defaultdict(set)
    indegree = [0] * len(sequence)
    for bond in mol.bonds:
        left, right = owner[bond.parent], owner[bond.child]
        if left == right:
            continue
        earlier, later = sorted((left, right))
        if later in successors[earlier]:
            continue
        successors[earlier].add(later)
        indegree[later] += 1

    order: list[int] = []
    available = [index for index, degree in enumerate(indegree) if degree == 0]
    while available:
        current = max(available)  # Prefer a different legal order when possible.
        available.remove(current)
        order.append(current)
        for successor in successors[current]:
            indegree[successor] -= 1
            if indegree[successor] == 0:
                available.append(successor)
    assert len(order) == len(sequence)
    return order


def _local_component(
    mol: Molecule,
    sequence: list[tuple[frozenset[int], ElemType]],
    block_index: int,
) -> Molecule:
    block = sequence[block_index][0]
    owner = {atom: index for index, (piece, _) in enumerate(sequence) for atom in piece}
    atoms = {
        atom: Atom(atom, mol.atoms[atom].layer, list(mol.atoms[atom].ends))
        for atom in block
    }
    bonds: list[Bond] = []
    for bond in mol.bonds:
        parent_in = bond.parent in block
        child_in = bond.child in block
        if parent_in and child_in:
            bonds.append(bond)
            continue
        if parent_in == child_in:
            continue

        inside = bond.parent if parent_in else bond.child
        outside = bond.child if parent_in else bond.parent
        direction = Direction.BOTTOM if parent_in else Direction.TOP
        kind = EndKind.FREE if block_index < owner[outside] else EndKind.FIXED
        atoms[inside].ends.append((kind, direction))

    component = Molecule(atoms, bonds)
    component.check_invariants()
    return component


def _assert_layout_replays(mol: Molecule) -> None:
    for sequence in _complete_sequences(mol):
        for index, (block, expected_kind) in enumerate(sequence):
            component = _local_component(mol, sequence, index)
            assert classify(component, set(block)) is expected_kind

        state = mol
        for index in _quotient_order(sequence, mol):
            block, expected_kind = sequence[index]
            assert classify(state, set(block)) is expected_kind
            state = state.cut_as_free(set(block))
        assert not state.atoms


def _assert_sequence_layout_replays(
    mol: Molecule, sequence: list[tuple[frozenset[int], ElemType]]
) -> None:
    for index, (block, expected_kind) in enumerate(sequence):
        component = _local_component(mol, sequence, index)
        assert classify(component, set(block)) is expected_kind

    state = mol
    for index in _quotient_order(sequence, mol):
        block, expected_kind = sequence[index]
        assert classify(state, set(block)) is expected_kind
        state = state.cut_as_free(set(block))
    assert not state.atoms


def _partitions(ids: frozenset[int]) -> list[tuple[frozenset[int], ...]]:
    """Enumerate set partitions with blocks of size at most three once each."""
    if not ids:
        return [()]

    first = min(ids)
    rest = ids - {first}
    out: list[tuple[frozenset[int], ...]] = []
    for size in range(1, min(3, len(ids)) + 1):
        for tail in combinations(sorted(rest), size - 1):
            block = frozenset({first, *tail})
            for suffix in _partitions(ids - block):
                out.append((block, *suffix))
    return out


def _layout_order(
    blocks: tuple[frozenset[int], ...], mol: Molecule, directions: dict[tuple[int, int], int]
) -> list[int] | None:
    owner = {atom: index for index, block in enumerate(blocks) for atom in block}
    successors: dict[int, set[int]] = defaultdict(set)
    indegree = [0] * len(blocks)
    for bond in mol.bonds:
        left, right = owner[bond.parent], owner[bond.child]
        if left == right:
            continue
        pair = tuple(sorted((left, right)))
        source = directions[pair]
        target = right if source == left else left
        if target not in successors[source]:
            successors[source].add(target)
            indegree[target] += 1

    order: list[int] = []
    available = [index for index, degree in enumerate(indegree) if degree == 0]
    while available:
        current = available.pop()
        order.append(current)
        for successor in successors[current]:
            indegree[successor] -= 1
            if indegree[successor] == 0:
                available.append(successor)
    return order if len(order) == len(blocks) else None


def _layout_component(
    mol: Molecule,
    blocks: tuple[frozenset[int], ...],
    block_index: int,
    directions: dict[tuple[int, int], int],
) -> Molecule:
    block = blocks[block_index]
    owner = {atom: index for index, piece in enumerate(blocks) for atom in piece}
    atoms = {
        atom: Atom(atom, mol.atoms[atom].layer, list(mol.atoms[atom].ends))
        for atom in block
    }
    bonds: list[Bond] = []
    for bond in mol.bonds:
        parent_in = bond.parent in block
        child_in = bond.child in block
        if parent_in and child_in:
            bonds.append(bond)
            continue
        if parent_in == child_in:
            continue

        inside = bond.parent if parent_in else bond.child
        outside = bond.child if parent_in else bond.parent
        pair = tuple(sorted((block_index, owner[outside])))
        direction = Direction.BOTTOM if parent_in else Direction.TOP
        kind = EndKind.FREE if directions[pair] == block_index else EndKind.FIXED
        atoms[inside].ends.append((kind, direction))

    component = Molecule(atoms, bonds)
    component.check_invariants()
    return component


def _assert_every_small_layout_replays(mol: Molecule) -> None:
    """Check the reverse G1 direction without deriving layouts from runs."""
    for blocks in _partitions(frozenset(mol.atoms)):
        owner = {atom: index for index, block in enumerate(blocks) for atom in block}
        pairs = sorted(
            {
                tuple(sorted((owner[bond.parent], owner[bond.child])))
                for bond in mol.bonds
                if owner[bond.parent] != owner[bond.child]
            }
        )
        for choices in product((0, 1), repeat=len(pairs)):
            directions = {
                pair: pair[choice]
                for pair, choice in zip(pairs, choices, strict=True)
            }
            order = _layout_order(blocks, mol, directions)
            if order is None:
                continue
            kinds = [
                classify(_layout_component(mol, blocks, index, directions), set(block))
                for index, block in enumerate(blocks)
            ]
            if ElemType.INVALID in kinds:
                continue

            state = mol
            for index in order:
                block = blocks[index]
                assert classify(state, set(block)) is kinds[index]
                state = state.cut_as_free(set(block))
            assert not state.atoms


@pytest.mark.parametrize("n", range(1, 4))
def test_layout_replays_every_small_toy1_sequence(n: int) -> None:
    for molecule in enumerate_toy1(n):
        _assert_layout_replays(molecule)


def test_layout_replays_a_branched_upper_layer() -> None:
    molecule = build(
        3,
        3,
        [(0, 2), (1, 2)],
        [(0, 1), (2, 1)],
        [(0, 0), (1, 0), (2, 2)],
    )
    _assert_layout_replays(molecule)


@pytest.mark.parametrize("n", range(1, 3))
def test_every_small_locally_legal_layout_replays(n: int) -> None:
    for molecule in enumerate_toy1(n):
        _assert_every_small_layout_replays(molecule)


@pytest.mark.parametrize("builder", [transpose_family, padded_family, trimmed_family])
def test_layout_replays_good_and_macro_block_shapes(builder) -> None:
    parameter = 2 if builder is trimmed_family else 1
    molecule = builder(parameter)
    record = run_toy1_algorithm(molecule)
    assert record.failed is None
    sequence = [(step.subset, step.elem_types[0]) for step in record.steps]

    kinds = {kind for _, kind in sequence}
    assert ElemType.E343 in kinds
    if builder is not transpose_family:
        assert ElemType.E33 in kinds
    _assert_sequence_layout_replays(molecule, sequence)
