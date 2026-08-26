# osmordred

**0D / 1D / 2D molecular descriptors in C++, on top of the RDKit.** Inspired by
Dragon, PaDEL and Mordred, reimplemented against the RDKit at the C++ level so a
full descriptor block is fast and parallelisable.

Around 3,500 individual features across 61 descriptor blocks, computed from a
`Chem.Mol` with no Python-side descriptor code in the hot path.

The exact width depends on the build and on the `version` argument, so check it
rather than hardcoding it — `len(calc_osmordred("CCCO"))` with the recipe below.
Earlier revisions of this file quoted 3586; the RDKit-integrated builds report
3585 (v2) and 3588 (v3). A feature-name list is only valid for the build that
produced it.

> **Read this first.** There is no `CalcOsmordred(mol)`. The package exposes one
> function *per descriptor block*, and the "compute everything" routine is a
> Python function you assemble yourself — the reference implementation is
> [`skbuild/test/test.py`](skbuild/test/test.py), reproduced under
> [Computing every descriptor](#computing-every-descriptor) below.

---

## Build

Two routes, for two different things. Both need Python 3.11 and RDKit 2023.9.3.

| Route | What you get | Where |
|---|---|---|
| **A — standalone wheel** | `import osmordred` as its own module, alongside a stock RDKit | `skbuild/` |
| **B — patched RDKit wheel** | Descriptors compiled *into* an RDKit wheel (`rdkit-pypi` + patches) | `build_patched_rdkit-pypi/` |

Route B is what the `cibuildwheel` GitHub workflow runs (manual dispatch,
macOS/`cp311-macosx_arm64`).

### Route A — standalone wheel

**All commands run from `skbuild/`, not the repository root.**

```bash
cd skbuild
./setup_env.sh          # creates conda env "osmordred" with the pinned deps
./build.sh              # conda run -n osmordred python -m build
conda activate osmordred
pip install dist/osmordred-0.2.0-*.whl --force-reinstall
```

`setup_env.sh` pins `boost==1.82.0 eigen lapack ninja python-build rdkit==2023.9.3`
on Python 3.11, and picks the BLAS for your platform — **openblas on macOS, MKL on
Linux**. It removes any existing `osmordred` env first.

To build into an environment you already have, do the same install by hand:

```bash
cd skbuild
rm -rf _skbuild dist
conda install -c conda-forge boost==1.82.0 eigen lapack ninja python-build \
    rdkit==2023.9.3 'blas=*=*openblas'      # Linux: 'blas=*=*mkl'
conda run python -m build
pip install dist/osmordred-0.2.0-*.whl --force-reinstall
pip show osmordred
```

### Route B — patched RDKit wheel

```bash
cd build_patched_rdkit-pypi
CIBW_PLATFORM=linux CIBW_BUILD=cp311-manylinux_x86_64 ./build_rdkit-pypi.sh
```

Wheels land in `wheelhouse/`. The patch set applied to `rdkit-pypi` 2023.09.3
lives in [`osmordred_rdkit_2023_09_3_patches/`](osmordred_rdkit_2023_09_3_patches/).

### Verifying the build

```bash
cd skbuild/test
pip install tqdm
python test.py          # writes Myfeatures.csv
```

**`import osmordred` fails from inside `skbuild/`** — Python picks up the
`osmordred/` *source subdirectory* instead of the installed package. Run from
anywhere else.

---

## Usage

Every descriptor block is its own function taking an RDKit `Mol` and returning a
sequence of doubles:

```python
from rdkit import Chem
import osmordred as rd

mol = Chem.MolFromSmiles("CCCO")

list(rd.CalcBalabanJ(mol))
list(rd.CalcInformationContent(mol, 5))     # depth argument
list(rd.CalcAtomCount(mol, 2))              # version argument — see below
```

### Computing every descriptor

Concatenate the blocks. This is the reference recipe from
[`skbuild/test/test.py`](skbuild/test/test.py), condensed:

```python
import numpy as np
from rdkit import Chem
import osmordred as rd

def calc_osmordred(smiles, version=2):
    """Every descriptor block, concatenated. version=1 reproduces Mordred's set."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    v = 1 if version == 1 else 2
    do_ex_estate = version != 1

    out = [
        rd.CalcABCIndex(mol), rd.CalcAcidBase(mol),
        rd.CalcAdjacencyMatrix(mol, v), rd.CalcAromatic(mol),
        rd.CalcAtomCount(mol, v), rd.CalcAutocorrelation(mol),
        # ... 55 more blocks, in the order test.py uses ...
        rd.CalcEState(mol, do_ex_estate),
        rd.CalcCarbonTypes(mol, v), rd.CalcDistanceMatrix(mol, v),
    ]
    if version != 1:
        out += [rd.CalcLogS(mol), rd.CalcInformationContent(mol, 5)]

    return np.concatenate([np.atleast_1d(np.array(x)) for x in out])
```

**Block order is part of the contract.** Concatenation is positional, so any
feature-name mapping you build is tied to the order you called the blocks in.
Take the order from `test.py` rather than from this README's abbreviated list.

### The `version` argument

Several blocks take a trailing `version` int, and the recipe threads one value
through all of them:

| | `version=1` | `version=2` (default) |
|---|---|---|
| Intent | Reproduce Mordred's original set | Extended set with the C++ fixes |
| `CalcAdjacencyMatrix` | — | adds `sm1` |
| `CalcAtomCount` | — | adds `nHetero` |
| `CalcCarbonTypes` | — | adds `FractionCSP3` |
| `CalcEState(mol, doExEstate)` | `False` | `True` |
| `CalcLogS` | not included | included |
| `CalcInformationContent(mol, 5)` | not included | included |

So **the two versions do not return the same width**, and a feature-name list
built against one will not resolve against the other.

---

## Descriptor blocks

The 61 blocks the reference recipe concatenates:

| Group | Blocks |
|---|---|
| Counts & constitution | `CalcAtomCount`, `CalcBondCount`, `CalcRingCount`, `CalcConstitutional`, `CalcCarbonTypes`, `CalcAromatic`, `CalcFrags`, `CalcFramework` |
| Physicochemical | `CalcWeight`, `CalcSLogP`, `CalcTopoPSA`, `CalcLipinski`, `CalcHydrogenBond`, `CalcAcidBase`, `CalcPolarizability`, `CalcMcGowanVolume`, `CalcVdwVolumeABC`, `CalcLogS`, `CalcMR`, `CalcPol`, `CalcODT`, `CalcAbrahams`, `CalcRotatableBond`, `CalcFlexibility` |
| Topological indices | `CalcABCIndex`, `CalcBalabanJ`, `CalcBertzCT`, `CalcWienerIndex`, `CalcZagrebIndex`, `CalcSchultz`, `CalcTopologicalIndex`, `CalcTopologicalCharge`, `CalcEccentricConnectivityIndex`, `CalcVertexAdjacencyInformation`, `CalcFragmentComplexity`, `CalcMolecularId` |
| Matrix-derived | `CalcAdjacencyMatrix`, `CalcDistanceMatrix`, `CalcDetourMatrix`, `CalcBaryszMatrix`, `CalcANMat`, `CalcASMat`, `CalcAZMat`, `CalcDSMat`, `CalcDN2Mat` |
| Connectivity & paths | `CalcChi`, `CalcPathCount`, `CalcWalkCount`, `CalcKappaShapeIndex`, `CalcAlphaKappaShapeIndex`, `CalcMolecularDistanceEdge` |
| Electrotopological | `CalcEState`, `CalcHEState`, `CalcBEState`, `CalcMoeType`, `CalcRNCGRPCG` |
| Autocorrelation & spectral | `CalcAutocorrelation`, `CalcBCUT` |
| Information content | `CalcInformationContent` |
| ETA | `CalcExtendedTopochemicalAtom` |
| Extras | `CalcAddFeatures` |

A further **34 functions are importable but not in the recipe** — the `*Eigen`
variants (`CalcAdjacencyMatrixEigen`, `CalcBaryszMatrixEigen`,
`CalcDetourMatrixEigen`, `CalcDistanceMatrixEigen`), the individual matrix
components (`CalcANI`, `CalcASN`, `CalcAZV`, `CalcDN2Z`, …), the Chi sub-blocks
(`CalcChiPath`, `CalcChichain`, `CalcChicluster`, `CalcChipathcluster`),
`CalcAcidicGroupCount`, `CalcBasicGroupCount`, `CalcCountAromaticAtoms`,
`CalcCountAromaticBonds` and `CalcExtendedTopochemicalAtom2`.
`dir(osmordred)` lists all 95.

Three more are defined in the C++ wrapper but **not re-exported** by
`osmordred/__init__.py`, so they are reachable only as
`osmordred._osmordred.<name>`: `CalcAutocorrelationEigen`, `CalcBCUTEigen` and
`CalcWalkCountBlas`. 98 defined, 95 importable.

---

## Accuracy and reproducibility

**Information Content is deliberately not bit-identical to Mordred or PaDEL.**
The implementation follows the original 1984 Basak method
([doi:10.1016/B978-0-08-030156-3.50138-7](https://doi.org/10.1016/B978-0-08-030156-3.50138-7))
within the RDKit, where aromaticity perception is its own thing — so it follows
Basak's core logic rather than reproducing another toolkit's output. The Basak
triplet features were added in the same effort.

**LAPACK, not Eigen3.** LAPACK was chosen for speed, particularly on the SVD of
symmetric square matrices. The cost is that it can produce very small
fluctuations on near-zero eigenvalues, which very slightly affects a few
descriptors. If you need bit-stable reruns, this is the reason you may not get
them.

**Parallelism** is handled in C++, so a batch does not pay the GIL.

## Known issues

This molecule causes a long stall — the computation is genuinely intensive
rather than hung:

```
c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41
```

The descriptor implementation is currently ~10k lines in a single file, and the
Python bindings would benefit from being split up and refactored.

## Repository layout

| Path | What it is |
|---|---|
| `Code/Osmordred/` | The C++ implementation (`Osmordred.cpp`/`.h`) and its Boost.Python wrapper (`Wrap/rdOsmordred.cpp`) |
| `skbuild/` | Route A — the standalone `osmordred` Python package, its build scripts and tests |
| `build_patched_rdkit-pypi/` | Route B — script and diff to build a patched `rdkit-pypi` wheel |
| `osmordred_rdkit_2023_09_3_patches/` | The patches Route B applies to RDKit 2023.09.3 |
| `rdkit/Chem/Osmordred/` | Python-side RDKit package files, for the Route B patched build |

## License

BSD-3-Clause
