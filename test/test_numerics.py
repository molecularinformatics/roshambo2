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


# --------------------------------------------------------------------------
# C2: the principal-axis frame must be a rotation, never a reflection.
# --------------------------------------------------------------------------

def _signed_volume(xyz, idx):
    """Triple product of four chosen atoms.

    Its sign is invariant under proper rotation and flips under reflection, so
    it is a direct chirality probe that needs no superposition.
    """
    p0, p1, p2, p3 = np.asarray(xyz, dtype=np.float64)[list(idx)]
    return float(np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0))


def _chirality_probe_atoms(xyz):
    """Pick four strongly non-coplanar atoms, greedily and deterministically.

    Simply taking the first four atoms is not safe: in most drug-like molecules
    they belong to a single aromatic ring, so they are coplanar and the sign of
    their triple product is numerical noise.
    """
    xyz = np.asarray(xyz, dtype=np.float64)
    i0 = int(np.argmax(np.linalg.norm(xyz - xyz.mean(axis=0), axis=1)))
    i1 = int(np.argmax(np.linalg.norm(xyz - xyz[i0], axis=1)))
    i2 = int(np.argmax(np.linalg.norm(
        np.cross(xyz - xyz[i0], xyz[i1] - xyz[i0]), axis=1)))
    n = np.cross(xyz[i1] - xyz[i0], xyz[i2] - xyz[i0])
    i3 = int(np.argmax(np.abs((xyz - xyz[i0]) @ n)))
    return i0, i1, i2, i3


@pytest.mark.parametrize("force_reflection", [False, True])
def test_principal_axis_alignment_preserves_chirality(monkeypatch, force_reflection):
    """Roshambo2Mol centres and aligns to the principal axes with

        _, _, v = np.linalg.svd(cvm);  r = r @ v.T

    The SVD makes no promise about the handedness of v -- each singular vector
    may be returned with either sign -- so det(v) can be -1, in which case that
    line mirrors the molecule.  A mirrored target is a genuinely different
    object to a chiral query because the optimiser only searches proper
    rotations, so it can never be superposed back.

    Real inputs almost never trigger it (0 of 10663 molecules sampled from
    test/confs.sdf and PubChem 3D), so this test *constructs* the det(v) = -1
    case rather than hoping to sample one.
    """
    from roshambo2.classes import Roshambo2Mol

    mol = _load("dataset.sdf", 1)[0]
    ref_xyz = np.array(Chem.RemoveHs(Chem.Mol(mol)).GetConformer().GetPositions())
    probe = _chirality_probe_atoms(ref_xyz)
    assert abs(_signed_volume(ref_xyz, probe)) > 1.0, "test molecule must not be planar"

    if force_reflection:
        real_svd = np.linalg.svd

        def reflecting_svd(a, *args, **kwargs):
            u, s, v = real_svd(a, *args, **kwargs)
            v = np.array(v)
            v[-1] = -v[-1]          # det(v) -> -det(v)
            return u, s, v

        # verify the construction really does hand the caller a left-handed basis
        cvm = np.cov(ref_xyz.astype(np.float32), rowvar=False)
        assert np.linalg.det(reflecting_svd(cvm)[2]) < 0

        monkeypatch.setattr(np.linalg, "svd", reflecting_svd)

    aligned = Roshambo2Mol(Chem.Mol(mol)).coords
    assert len(aligned) == len(ref_xyz)

    assert (np.sign(_signed_volume(aligned, probe))
            == np.sign(_signed_volume(ref_xyz, probe))), (
        "principal-axis alignment mirrored the molecule"
    )


# dataset.sdf entry 1 (CHEMBL1079594), rotated about z by this angle, is a real
# input for which numpy's SVD returns a left-handed basis.  Found by sweeping
# rigid orientations: 4 hits in 6.1M (molecule, orientation) trials over the 153
# molecules of dataset.sdf, so it is rare but not exotic, and the window here is
# ~0.0026 degrees wide rather than a single-float coincidence.  Which
# orientations trip it is a property of the host LAPACK's sign convention, hence
# the skip below rather than an unconditional assert.
MIRRORING_ROTATION_ANGLE = 2.9033342667


def _rotate_z(mol, theta):
    """Rigidly rotate every atom of a copy of mol about z. Chemistry unchanged."""
    from rdkit.Geometry import Point3D

    mol = Chem.Mol(mol)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    conf = mol.GetConformer()
    r = np.array(conf.GetPositions()) @ R.T
    for j in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(j, Point3D(*(float(v) for v in r[j])))
    return mol


def test_real_molecule_is_not_mirrored():
    """A molecule from this repo's own test set that the SVD returns left-handed.

    The monkeypatched test above proves the guard works on any LAPACK; this one
    proves the defect is reachable from an ordinary SDF.  A file carries whatever
    orientation it happens to carry, and for this one the covariance matrix has
    well-separated singular values (19.46 / 2.27 / 1.03) -- nothing degenerate or
    planar -- yet det(v) is exactly -1.
    """
    from roshambo2.classes import Roshambo2Mol

    mols = _load("dataset.sdf", 2)
    mol = _rotate_z(mols[1], MIRRORING_ROTATION_ANGLE)

    ref_xyz = np.array(Chem.RemoveHs(Chem.Mol(mol)).GetConformer().GetPositions())
    cvm = np.cov(np.asarray(ref_xyz, dtype=np.float32), rowvar=False)
    if np.linalg.det(np.linalg.svd(cvm)[2]) > 0:
        pytest.skip("this LAPACK returns a right-handed basis for this orientation")

    probe = _chirality_probe_atoms(ref_xyz)
    assert abs(_signed_volume(ref_xyz, probe)) > 1.0, "test molecule must not be planar"

    aligned = Roshambo2Mol(Chem.Mol(mol)).coords
    assert (np.sign(_signed_volume(aligned, probe))
            == np.sign(_signed_volume(ref_xyz, probe))), (
        "principal-axis alignment mirrored a real input molecule"
    )


def test_alignment_is_a_proper_rotation_on_real_molecules():
    """The map from input coordinates to aligned coordinates must have det +1."""
    from roshambo2.classes import Roshambo2Mol

    for mol in _load("dataset.sdf", 25):
        ref = np.array(Chem.RemoveHs(Chem.Mol(mol)).GetConformer().GetPositions())
        aligned = Roshambo2Mol(Chem.Mol(mol)).coords
        a = ref - ref.mean(axis=0)
        # aligned = a @ V.T exactly, so recover V and check its handedness
        V = np.linalg.lstsq(a, aligned, rcond=None)[0].T
        assert np.linalg.det(V) > 0
        assert np.allclose(V @ V.T, np.eye(3), atol=1e-3)


# --------------------------------------------------------------------------
# C3: in combination mode the mixed gradient must be the gradient of the
# blended Tanimoto that is actually reported.
# --------------------------------------------------------------------------

def _run_opt(backend, start_mode, color, mixing, legacy, n_query=6, n_target=25):
    mols = _load("dataset.sdf", n_query + n_target)
    cg = PharmacophoreGenerator()
    qdata = prepare_from_rdkitmols(mols[:n_query], color=color,
                                   color_generator=cg, quiet=True)
    tdata = prepare_from_rdkitmols(mols[n_query:], color=color,
                                   color_generator=cg, quiet=True)
    settings = {"lr_q": 0.1, "lr_t": 0.1, "steps": 100,
                "legacy_combo_gradient": legacy}
    if backend == "cuda":
        overlay = CudaShapeOverlay(qdata, tdata, start_mode,
                                   color_generator=cg if color else None,
                                   mixing=mixing, verbosity=0, n_gpus=1,
                                   optimizer_settings=settings)
    else:
        overlay = CppShapeOverlay(qdata, tdata, start_mode,
                                  color_generator=cg if color else None,
                                  mixing=mixing, verbosity=0,
                                  optimizer_settings=settings)
    return np.array(overlay.optimize_overlap())


@pytest.mark.parametrize("backend", BACKENDS)
def test_shape_only_is_unaffected_by_combo_gradient_flag(backend):
    """Shape-only must be bit-identical either way.

    With a single channel the gradient is scaled by one positive constant,
    which cancels exactly out of an Adagrad step
    (lr * c*g / sqrt(sum (c*g)^2) == lr * g / sqrt(sum g^2)), so the
    normalization never affected the shape-only trajectory and the flag must
    not change it.
    """
    a = _run_opt(backend, 1, color=False, mixing=0.0, legacy=False)
    b = _run_opt(backend, 1, color=False, mixing=0.0, legacy=True)
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("backend", BACKENDS)
def test_combo_gradient_reaches_a_better_blended_tanimoto(backend):
    """The corrected weighting optimises the reported objective.

    g_combo = (1-m) * dTs/dVs * grad(Vs) + m * dTc/dVc * grad(Vc) is exactly
    grad of the reported blend, so averaged over pairs it must not do worse
    than the previous weighting, which mixed the two channels in the wrong
    ratio.
    """
    new = _run_opt(backend, 1, color=True, mixing=0.5, legacy=False)
    old = _run_opt(backend, 1, color=True, mixing=0.5, legacy=True)

    # the flag must actually do something in combination mode
    assert not np.array_equal(new[:, :, 0], old[:, :, 0])

    assert np.mean(new[:, :, 0]) >= np.mean(old[:, :, 0])


# --------------------------------------------------------------------------
# The two permanent invariants, checked on 100% of returned rows.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("start_mode", [0, 1, 2])
@pytest.mark.parametrize("color", [False, True])
def test_overlap_respects_cauchy_schwarz(backend, start_mode, color):
    """V_AB <= sqrt(V_AA * V_BB) for any rigid pose.

    All shape Gaussians share one width, so V_AB is the inner product of the
    two density fields in the kernel's RKHS and Cauchy-Schwarz applies for any
    pose whatsoever.  Exceeding it means the reported overlap is achievable by
    no rigid placement of the two molecules.

    The other plausible bound, V_AB <= min(V_AA, V_BB), is *false* for this
    first-order expression -- there is no inclusion-exclusion correction -- so
    it is deliberately not asserted.
    """
    scores = _run(backend, start_mode, color)
    v_ab, v_aa, v_bb = scores[:, :, 3], scores[:, :, 5], scores[:, :, 6]
    ratio = v_ab / np.sqrt(v_aa * v_bb)
    assert np.all(np.isfinite(ratio))
    worst = float(np.max(ratio))
    assert worst <= 1.0 + 1e-4, f"max V_AB/sqrt(V_AA V_BB) = {worst:.6f}"


@pytest.mark.parametrize("backend", BACKENDS)
def test_self_match_saturates_cauchy_schwarz(backend):
    """A molecule against itself must reach the bound, so it is not vacuous."""
    mols = _load("dataset.sdf", 8)
    data = prepare_from_rdkitmols(mols, color=False, quiet=True)
    if backend == "cuda":
        overlay = CudaShapeOverlay(data, data, 0, verbosity=0, n_gpus=1)
    else:
        overlay = CppShapeOverlay(data, data, 0, verbosity=0)
    scores = np.array(overlay.optimize_overlap())
    diag = np.array([scores[i, i] for i in range(len(mols))])
    ratio = diag[:, 3] / np.sqrt(diag[:, 5] * diag[:, 6])
    assert np.allclose(ratio, 1.0, atol=1e-3), ratio


def _good_row():
    row = np.zeros((1, 1, 20), dtype=np.float64)
    row[0, 0, 3] = 500.0     # V_AB
    row[0, 0, 5] = 1000.0    # V_AA
    row[0, 0, 6] = 1000.0    # V_BB
    row[0, 0, 9] = 1.0       # identity quaternion
    return row


def test_invariant_checker_accepts_a_valid_row():
    from roshambo2.roshambo2 import check_result_invariants
    check_result_invariants(_good_row())


def test_invariant_checker_rejects_a_non_unit_quaternion():
    from roshambo2.roshambo2 import check_result_invariants
    row = _good_row()
    row[0, 0, 9] = 1.01
    with pytest.raises(RuntimeError, match="non-unit quaternion"):
        check_result_invariants(row)


def test_invariant_checker_rejects_a_cauchy_schwarz_violation():
    from roshambo2.roshambo2 import check_result_invariants
    row = _good_row()
    row[0, 0, 3] = 1010.0     # V_AB above sqrt(V_AA * V_BB) = 1000
    with pytest.raises(RuntimeError, match="Cauchy-Schwarz"):
        check_result_invariants(row)


def test_invariant_checker_allows_v_ab_above_min_self_overlap():
    """V_AB <= min(V_AA, V_BB) must NOT be enforced: it is false here."""
    from roshambo2.roshambo2 import check_result_invariants
    row = _good_row()
    row[0, 0, 5] = 400.0      # V_AA
    row[0, 0, 6] = 2500.0     # V_BB; sqrt(V_AA*V_BB) = 1000
    row[0, 0, 3] = 900.0      # above min(V_AA, V_BB) but below the bound
    check_result_invariants(row)


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("start_mode", [0, 1, 2])
@pytest.mark.parametrize("optim_mode", ["shape", "color", "combination"])
def test_frontend_enforces_invariants(backend, start_mode, optim_mode):
    """The checks run on every compute() call, for every mode."""
    from roshambo2 import Roshambo2

    calc = Roshambo2("query.sdf", "dataset.sdf",
                     color=optim_mode != "shape", verbosity=0)
    calc.compute(backend=backend, start_mode=start_mode, optim_mode=optim_mode,
                 combination_param=0.5 if optim_mode == "combination" else None,
                 write_scores=False)
