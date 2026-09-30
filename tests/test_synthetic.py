"""Independent algebraic checks of the frozen synthetic constructions."""

import json

import anndata as ad
import numpy as np
import pytest

from vec_scorelens.synthetic import CASES, create_case, create_healthy_sample, write_case


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("name", CASES)
def test_fixtures_are_finite_nonnegative_reproducible_and_unlabeled(name, seed):
    case = create_case(name, seed)
    assert case["planted"]["synthetic_only"] is True
    assert case["planted"]["matrix_sha256"] == create_case(name, seed)["planted"]["matrix_sha256"]
    for key in ("prediction", "target", "reference"):
        assert np.isfinite(case[key].X).all()
        assert np.all(case[key].X >= 0)
        assert case[key].n_vars == 128
        assert case[key].var_names.is_unique
    assert "celltype" in case["target"].obs
    assert "celltype" not in case["prediction"].obs
    assert "celltype" not in case["reference"].obs


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_exact_response_and_sign_magnitude_constructions(seed):
    for name, expected in (("direction", -1.0), ("magnitude", 1.8)):
        c = create_case(name, seed)
        t, r, p = (c[k].X for k in ("target", "reference", "prediction"))
        dt, dp = t.mean(0) - r.mean(0), p.mean(0) - r.mean(0)
        selected = [64, 65, 66, 67, 80, 81, 82, 83]
        np.testing.assert_allclose(dp[selected], expected * dt[selected], atol=1e-10)
        other = np.setdiff1d(np.arange(128), selected)
        np.testing.assert_allclose(dp[other], dt[other], atol=1e-10)
    c = create_case("response_magnitude", seed)
    dt = c["target"].X.mean(0) - c["reference"].X.mean(0)
    dp = c["prediction"].X.mean(0) - c["reference"].X.mean(0)
    np.testing.assert_allclose(dp, 2 * dt, atol=1e-10)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_covariance_corruption_preserves_marginals_and_destroys_joint(seed):
    c = create_case("covariance", seed)
    p, t = c["prediction"].X, c["target"].X
    np.testing.assert_allclose(p.mean(0), t.mean(0), atol=1e-10)
    np.testing.assert_allclose(p.var(0), t.var(0), atol=1e-10)
    np.testing.assert_array_equal(np.sort(p, axis=0), np.sort(t, axis=0))
    truth_correlations, predicted_correlations = [], []
    for module in c["planted"]["covariance_modules"]:
        a = np.corrcoef(t[:, module], rowvar=False)
        b = np.corrcoef(p[:, module], rowvar=False)
        idx = np.triu_indices(8, 1)
        truth_correlations.extend(np.abs(a[idx]))
        predicted_correlations.extend(np.abs(b[idx]))
    assert np.median(truth_correlations) >= 0.80
    assert np.median(predicted_correlations) <= 0.20


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_population_and_state_shift_ground_truth(seed):
    c = create_case("missing", seed)
    assert c["prediction"].n_obs == 864
    assert "state_3" not in c["planted"]["privileged_test_only"]["prediction_state_ids"]
    c = create_case("overrepresented", seed)
    labels = np.asarray(c["planted"]["privileged_test_only"]["prediction_state_ids"])
    assert np.sum(labels == "state_3") == 288
    c = create_case("within_state", seed)
    labels = np.asarray(c["planted"]["privileged_test_only"]["prediction_state_ids"])
    difference = c["prediction"].X - c["target"].X
    expected = np.zeros_like(difference)
    expected[np.ix_(np.flatnonzero(labels == "state_2"), np.arange(96, 104))] = 0.75
    np.testing.assert_allclose(difference, expected, atol=1e-10)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_scale_collapse_and_healthy_invariants(seed):
    c = create_case("scale", seed)
    np.testing.assert_allclose(c["prediction"].X, 1.8 * c["target"].X + 0.7, atol=1e-10)
    c = create_case("collapse", seed)
    np.testing.assert_allclose(c["prediction"].X.mean(0), c["target"].X.mean(0), atol=1e-10)
    assert c["prediction"].X.var(0).mean() <= 1e-20
    c = create_healthy_sample(seed, 0)
    assert not np.array_equal(c["prediction"].X, c["target"].X)
    np.testing.assert_allclose(c["prediction"].X.mean(0), c["target"].X.mean(0), atol=1e-10)


def test_write_case_creates_named_inputs_and_refuses_existing_output(tmp_path):
    out = write_case(tmp_path / "new_demo", "direction", 1)
    assert {p.name for p in out.iterdir()} == {"prediction.h5ad", "target.h5ad", "reference.h5ad", "truth.json"}
    assert json.loads((out / "truth.json").read_text())["seed"] == 1
    for name in ("prediction", "target", "reference"):
        assert ad.read_h5ad(out / f"{name}.h5ad").n_vars == 128
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    with pytest.raises(FileExistsError):
        write_case(out, "collapse", 2)
    assert before == {p.name: p.read_bytes() for p in out.iterdir()}


def test_unknown_case_and_negative_seed_are_rejected():
    with pytest.raises(ValueError):
        create_case("not_a_case")
    with pytest.raises(ValueError):
        create_case("healthy", -1)
