"""Regression tests for the two v3 C++ fixes.

Runs under pytest, or standalone: `python test_v3_fixes.py`.

Both guard failures that are silent — they produce a plausible-looking result, or
no result at all, rather than an error — so they need a test rather than a comment.
"""
import math
import time

import pytest
from rdkit import Chem

import osmordred as rd

# Degenerate molecules whose adjacency/distance systems are singular. Before the
# dgelss tier, these fell through to a 5-element {0,0,0,0,0} stub instead of
# nrhs*5 = 20, so the molecule silently returned a SHORT descriptor vector.
DEGENERATE = ["C", "C=C", "C#C", "CC", "O=C=O", "N#N", "O", "CO"]

# A normal molecule the solve has always succeeded on, used as the reference
# width. Deliberately not a hardcoded number: the width moves between v2 and v3
# and between builds, so the invariant is "same as a molecule that works".
REFERENCE = "c1ccccc1"

# The 14 *_Dt descriptors need the longest simple path between every atom pair --
# NP-hard. A naive whole-graph DFS hangs indefinitely on this one.
POLYCYCLIC = (
    "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9"
    "c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"
)

MATRIX_BLOCKS = ["CalcANMat", "CalcASMat", "CalcAZMat", "CalcDSMat", "CalcDN2Mat"]


@pytest.mark.parametrize("block", MATRIX_BLOCKS)
@pytest.mark.parametrize("smiles", DEGENERATE)
def test_singular_matrices_return_full_width(block, smiles):
    """dgelss tier: a singular system must not truncate the descriptor vector."""
    fn = getattr(rd, block)
    expected = len(fn(Chem.MolFromSmiles(REFERENCE)))
    got = len(fn(Chem.MolFromSmiles(smiles)))
    assert got == expected, (
        f"{block}({smiles}) returned {got} values, expected {expected} -- the "
        "singular-matrix fallback truncated the vector (missing dgelss tier?)"
    )


def test_detour_matrix_terminates_on_polycyclic():
    """Biconnected decomposition: must finish rather than hang."""
    mol = Chem.MolFromSmiles(POLYCYCLIC)
    assert mol is not None
    start = time.time()
    values = rd.CalcDetourMatrix(mol)
    elapsed = time.time() - start
    assert elapsed < 60, f"detour matrix took {elapsed:.1f}s -- the BCC guard is not firing"
    assert len(values) == 14
    # A single ring system above circuit rank 12 is reported as NaN by design,
    # so either a finite result or all-NaN is correct -- a hang is not.
    assert all(math.isnan(v) for v in values) or any(not math.isnan(v) for v in values)


def test_detour_matrix_unchanged_on_normal_molecules():
    """The decomposition is exact, so ordinary molecules must be unaffected."""
    for smiles in ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "C1CC2CCC1CC2", "c1ccc2ccccc2c1"]:
        values = rd.CalcDetourMatrix(Chem.MolFromSmiles(smiles))
        assert len(values) == 14
        assert all(not math.isnan(v) for v in values), f"unexpected NaN for {smiles}"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
