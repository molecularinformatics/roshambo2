"""Numerical-contract tests for the optimiser output.

These are cross-backend tests of properties that must hold for *every* returned
row, independent of which molecules were searched.  They are deliberately
stated as properties of the result rather than as comparisons between backends:
the pre-existing backend-vs-backend tests only check that the top-3 names
agree, which cannot catch a defect that both backends share.
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

QNORM_TOL = 1e-5


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


def _run(backend, start_mode, color, mixing=0.5, n_query=6, n_target=25):
    mols = _load("dataset.sdf", n_query + n_target)
    cg = PharmacophoreGenerator()
    qdata = prepare_from_rdkitmols(mols[:n_query], color=color,
                                   color_generator=cg, quiet=True)
    tdata = prepare_from_rdkitmols(mols[n_query:], color=color,
                                   color_generator=cg, quiet=True)
    if backend == "cuda":
        overlay = CudaShapeOverlay(qdata, tdata, start_mode,
                                   color_generator=cg if color else None,
                                   mixing=mixing if color else 0.0,
                                   verbosity=0, n_gpus=1)
    else:
        overlay = CppShapeOverlay(qdata, tdata, start_mode,
                                  color_generator=cg if color else None,
                                  mixing=mixing if color else 0.0,
                                  verbosity=0)
    return np.array(overlay.optimize_overlap())


def _quat_to_matrix(q):
    """The unit-normalised form the backends use to build the pose matrix."""
    w, x, y, z = q
    return np.array([
        [1 - 2 * y * y - 2 * z * z, 2 * x * y - 2 * w * z, 2 * x * z + 2 * w * y],
        [2 * x * y + 2 * w * z, 1 - 2 * x * x - 2 * z * z, 2 * y * z - 2 * w * x],
        [2 * x * z - 2 * w * y, 2 * y * z + 2 * w * x, 1 - 2 * x * x - 2 * y * y],
    ])


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("start_mode", [0, 1, 2])
@pytest.mark.parametrize("color", [False, True])
def test_returned_quaternion_is_unit(backend, start_mode, color):
    """The optimiser renormalises q every step; the result must be unit."""
    scores = _run(backend, start_mode, color)
    norms = np.linalg.norm(scores[:, :, 9:13], axis=-1)
    worst = float(np.max(np.abs(norms - 1.0)))
    assert worst <= QNORM_TOL, f"max | ||q|| - 1 | = {worst:.3e}"

    # the start-mode quaternion is built from an axis/angle and must also be unit
    start_norms = np.linalg.norm(scores[:, :, 16:20], axis=-1)
    assert np.max(np.abs(start_norms - 1.0)) <= QNORM_TOL


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("start_mode", [0, 1, 2])
def test_returned_pose_is_a_rotation(backend, start_mode):
    """quaternion_to_rotation_matrix assumes ||q|| == 1.

    Given a quaternion of true norm s it returns s**2 * R(q/s) + (1 - s**2) * I,
    which is a non-orthogonal anisotropic stretch rather than a rotation.  The
    pose the caller gets back -- and which get_best_fit_structures applies to
    the output coordinates -- must be a genuine rigid rotation: orthogonal,
    right-handed, unit determinant.
    """
    scores = _run(backend, start_mode, color=False)
    worst_orth = 0.0
    worst_det = 0.0
    for row in scores.reshape(-1, 20):
        M = _quat_to_matrix(row[9:13].astype(np.float64))
        worst_orth = max(worst_orth, float(np.max(np.abs(M.T @ M - np.eye(3)))))
        worst_det = max(worst_det, abs(float(np.linalg.det(M)) - 1.0))
    assert worst_orth < 1e-4, f"max ||M'M - I||_inf = {worst_orth:.3e}"
    assert worst_det < 1e-4, f"max |det(M) - 1| = {worst_det:.3e}"
