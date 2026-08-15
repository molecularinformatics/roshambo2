"""The colour output columns must be zero, not NaN, when colour scoring is off.

When colour is disabled the backends are handed a 1x1 all-zero interaction
matrix, but the colour loops still ran over whatever sits past the real atom
count.  Every such pair evaluated exp(-a*a*d2/(a+a)) and (pi/(2a))**1.5 with
a == 0, giving 0/0 and inf, so the colour volume came out NaN.
"""

import numpy as np
import pytest
from rdkit import Chem

from roshambo2.backends.cpp_backend import CppShapeOverlay
from roshambo2.pharmacophore import PharmacophoreGenerator
from roshambo2.prepare import prepare_from_rdkitmols

try:
    from roshambo2.backends.cuda_backend import CudaShapeOverlay
    HAVE_CUDA = True
except Exception:  # pragma: no cover - depends on the build
    HAVE_CUDA = False

BACKENDS = ["cpp"] + (["cuda"] if HAVE_CUDA else [])


def _load(path, n):
    mols = []
    for m in Chem.SDMolSupplier(path, removeHs=False):
        if m is None:
            continue
        m = Chem.Mol(m)
        m.SetProp("_Name", f"mol{len(mols)}")
        mols.append(m)
        if len(mols) >= n:
            break
    return mols


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("start_mode", [0, 1, 2])
@pytest.mark.parametrize("prepared_with_color", [False, True])
def test_shape_only_run_has_no_nans(backend, start_mode, prepared_with_color):
    """With colour scoring off the colour columns must be zero, not NaN.

    Two ways in, both covered here:
      * data prepared without colour -- the range past the real atom count holds
        the padding atoms, so any molecule set with ragged atom counts poisons
        overlap_color and self_overlap_color_* in the user-facing scores;
      * data prepared with colour but scored shape-only -- the range holds real
        colour features whose type indices also read past the end of the 1x1
        placeholder matrix.
    """
    mols = _load("dataset.sdf", 31)
    cg = PharmacophoreGenerator()
    kw = {"color_generator": cg} if prepared_with_color else {}
    qdata = prepare_from_rdkitmols(mols[:6], color=prepared_with_color, quiet=True, **kw)
    tdata = prepare_from_rdkitmols(mols[6:], color=prepared_with_color, quiet=True, **kw)

    if backend == "cuda":
        overlay = CudaShapeOverlay(qdata, tdata, start_mode, verbosity=0, n_gpus=1)
    else:
        overlay = CppShapeOverlay(qdata, tdata, start_mode, verbosity=0)
    scores = np.array(overlay.optimize_overlap())

    bad = sorted(set(np.argwhere(~np.isfinite(scores))[:, 2].tolist()))
    assert not bad, f"non-finite output columns: {bad}"
    assert np.all(scores[:, :, [4, 7, 8]] == 0.0)
