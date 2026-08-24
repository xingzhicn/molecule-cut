"""Finite regression guards for the G0 feedback witnesses."""

from __future__ import annotations

import pytest

from molecule_cut.enumerate import enumerate_toy1
from molecule_cut.gadgets import two_sets_yes
from molecule_cut.molecule import Bond, Molecule


def _is_forest(vertices: set[int], bonds: list[Bond]) -> bool:
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for bond in bonds:
        if bond.parent not in vertices or bond.child not in vertices:
            continue
        left, right = find(bond.parent), find(bond.child)
        if left == right:
            return False
        parent[left] = right
    return True


def _assert_g0_witnesses(mol: Molecule) -> None:
    cross = mol.cross_bonds()
    assert mol.circuit_rank() == max(len(cross) - 1, 0)

    # Keep one cross edge; all remaining cross edges form a feedback-edge set.
    removed_edges = set(cross[1:])
    retained_bonds = [bond for bond in mol.bonds if bond not in removed_edges]
    assert _is_forest(set(mol.atoms), retained_bonds)

    # One endpoint per removed edge is a feedback-vertex witness.
    feedback_vertices = {bond.parent for bond in removed_edges}
    assert len(feedback_vertices) <= mol.circuit_rank()
    assert _is_forest(set(mol.atoms) - feedback_vertices, mol.bonds)


@pytest.mark.parametrize("n", range(1, 5))
def test_g0_witnesses_on_toy1_enumeration(n: int) -> None:
    for molecule in enumerate_toy1(n):
        _assert_g0_witnesses(molecule)


def test_g0_witnesses_on_multicross_gadget() -> None:
    _assert_g0_witnesses(two_sets_yes())
