import pytest
from fastapi.testclient import TestClient

from app.config import settings
from app.enrichment import normalise_base_peak
from app.envelope import merge_peaks


@pytest.fixture()
def client(tmp_path, monkeypatch):
    # Isolate each test run with its own SQLite file.
    monkeypatch.setattr(settings, "db_path", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def _ion_peaks(client, formula, adduct, charge, top_n=30):
    resp = client.post(
        "/api/v1/ion-envelope",
        json={
            "formula": formula,
            "adduct": adduct,
            "charge": charge,
            "top_n": top_n,
        },
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    return [
        {"mz": p["mz"], "intensity": p["intensity"]} for p in body["peaks"]
    ]


def _mixed_peaks(peaks_a, peaks_b, fraction_a_percent, tolerance=0.01):
    """Reference implementation of the API's mixture construction."""
    fraction_a = fraction_a_percent / 100.0
    combined = [(p["mz"], p["intensity"] * fraction_a) for p in peaks_a]
    combined.extend(
        (p["mz"], p["intensity"] * (1.0 - fraction_a)) for p in peaks_b
    )
    merged = [(mz, i) for mz, i in merge_peaks(combined, tolerance) if i > 0]
    return [
        {"mz": round(mz, 6), "intensity": round(i, 6)}
        for mz, i in normalise_base_peak(merged)
    ]


def _payload(**overrides):
    payload = {
        "formula": "C6H12O6",
        "peaks": [{"mz": 181.07, "intensity": 100.0}],
        "mass_tolerance": 0.01,
        "candidates": [
            {"adduct": "+H", "charge": 1},
            {"adduct": "+Na-H", "charge": 1},
        ],
        "step_percent": 1.0,
    }
    payload.update(overrides)
    return payload


def test_ion_mixture_fit_picks_pair_and_ratio(client: TestClient):
    peaks_h = _ion_peaks(client, "C6H12O6", "+H", 1)
    peaks_na = _ion_peaks(client, "C6H12O6", "+Na-H", 1)
    observed = _mixed_peaks(peaks_h, peaks_na, 50.0)

    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(
            peaks=observed,
            step_percent=5.0,
            candidates=[
                {"adduct": "+H", "charge": 1},
                {"adduct": "+Na-H", "charge": 1},
                {"adduct": "+2H", "charge": 2},
            ],
        ),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()

    assert body["formula"] == "C6H12O6"
    assert body["mass_tolerance"] == 0.01
    assert body["step_percent"] == 5.0
    assert body["observed_peak_count"] == len(observed)
    assert body["valid_candidate_count"] == 3
    # 3 valid candidates -> 3 unordered pairs; 5% grid over 0-100% -> 21 points.
    assert body["combination_count"] == 3
    assert body["grid_points"] == 21
    assert body["search_points"] == 63

    # Combinations are ranked by ascending MSE.
    mses = [c["mse"] for c in body["combinations"]]
    assert mses == sorted(mses)
    assert [c["rank"] for c in body["combinations"]] == [1, 2, 3]

    best = body["best_mixture"]
    assert best["first_index"] == 0
    assert best["second_index"] == 1
    assert best["first_adduct"] == "+H"
    assert best["first_charge"] == 1
    assert best["second_adduct"] == "+Na-H"
    assert best["second_charge"] == 1
    assert best["first_fraction_percent"] == 50.0
    assert best["mse"] == pytest.approx(0.0, abs=0.05)

    winner = body["combinations"][0]
    assert winner["pair_count"] == len(winner["pairs"])
    # Every observed peak participates in exactly one serialised pair.
    observed_mzs = [
        p["observed_mz"] for p in winner["pairs"] if p["observed_mz"] is not None
    ]
    assert len(observed_mzs) == len(observed)
    assert set(winner["pairs"][0]) == {
        "observed_mz",
        "observed_intensity",
        "predicted_mz",
        "predicted_intensity",
    }
    # The full (non-truncated) predicted envelopes leave some tiny unmatched
    # predicted peaks with the observed side nulled out.
    unmatched = [p for p in winner["pairs"] if p["observed_mz"] is None]
    assert unmatched
    assert unmatched[0]["predicted_intensity"] is not None


def test_ion_mixture_fit_tie_breaks_by_index_pair_then_ratio(client: TestClient):
    # Three identical specifications give identical envelopes: every pair and
    # every mixing ratio produce the same MSE. The index-pair tie-break must
    # pick (0, 1); the ratio tie-break must pick the lowest first fraction.
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(
            step_percent=10.0,
            candidates=[
                {"adduct": "+H", "charge": 1},
                {"adduct": "+H", "charge": 1},
                {"adduct": "+H", "charge": 1},
            ],
        ),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["grid_points"] == 11
    pairs = [(c["first_index"], c["second_index"]) for c in body["combinations"]]
    assert pairs == [(0, 1), (0, 2), (1, 2)]
    best = body["best_mixture"]
    assert (best["first_index"], best["second_index"]) == (0, 1)
    assert best["first_fraction_percent"] == 0.0


def test_ion_mixture_fit_pure_endpoints_match_single_envelopes(client: TestClient):
    # A 0% first-candidate mixture must equal the second envelope; 100% the
    # first. Observed spectrum is the pure +Na-H envelope.
    peaks_na = _ion_peaks(client, "C6H12O6", "+Na-H", 1)
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(peaks=peaks_na, step_percent=10.0),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    best = body["best_mixture"]
    assert best["first_index"] == 0
    assert best["second_index"] == 1
    assert best["first_fraction_percent"] == 0.0
    assert best["mse"] == pytest.approx(0.0, abs=0.05)


def test_ion_mixture_fit_invalid_candidate_reports_and_skips_pairs(
    client: TestClient,
):
    peaks_h = _ion_peaks(client, "C6H12O6", "+H", 1)
    peaks_na = _ion_peaks(client, "C6H12O6", "+Na-H", 1)
    observed = _mixed_peaks(peaks_h, peaks_na, 30.0)

    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(
            peaks=observed,
            step_percent=10.0,
            candidates=[
                {"adduct": "+H", "charge": 1},        # valid, index 0
                {"adduct": "H", "charge": 1},         # missing sign, index 1
                {"adduct": "-99H", "charge": -1},     # negative count, index 2
                {"adduct": "+Na-H", "charge": 1},     # valid, index 3
            ],
        ),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()

    assert body["valid_candidate_count"] == 2
    # Only the pair of the two valid candidates (indices 0 and 3) is scored.
    assert body["combination_count"] == 1
    best = body["best_mixture"]
    assert (best["first_index"], best["second_index"]) == (0, 3)
    assert best["first_fraction_percent"] == 30.0
    assert best["mse"] == pytest.approx(0.0, abs=0.05)

    by_index = {c["index"]: c for c in body["candidates"]}
    assert [c["index"] for c in body["candidates"]] == [0, 1, 2, 3]
    assert by_index[0]["valid"] is True
    bad = by_index[1]
    assert bad["valid"] is False
    assert bad["error"] == "invalid_adduct"
    assert bad["position"] == 0
    bad = by_index[2]
    assert bad["valid"] is False
    assert bad["error"] == "invalid_adduct"
    assert bad["position"] == 3
    assert "negative" in bad["message"]
    assert by_index[3]["valid"] is True


def test_ion_mixture_fit_fewer_than_two_valid_candidates(client: TestClient):
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(
            candidates=[
                {"adduct": "+H", "charge": 1},
                {"adduct": "H", "charge": 1},
                {"adduct": "-99H", "charge": -1},
            ]
        ),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["valid_candidate_count"] == 1
    assert body["combination_count"] == 0
    assert body["grid_points"] == 0
    assert body["search_points"] == 0
    assert body["combinations"] == []
    assert body["best_mixture"] is None
    assert len(body["candidates"]) == 3
    assert [c["valid"] for c in body["candidates"]] == [True, False, False]


def test_ion_mixture_fit_unmatched_observed_peak_counts(client: TestClient):
    peaks_h = _ion_peaks(client, "C6H12O6", "+H", 1, top_n=3)
    observed = peaks_h + [
        {"mz": peaks_h[0]["mz"] + 50.0, "intensity": 40.0}
    ]
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(peaks=observed, step_percent=10.0),
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    best = body["best_mixture"]
    # The extra observed peak is far outside tolerance: its 40% normalised
    # intensity must make the best MSE clearly positive.
    assert best["mse"] > 1.0
    winner = body["combinations"][0]
    unmatched = [p for p in winner["pairs"] if p["predicted_mz"] is None]
    assert len(unmatched) == 1
    assert unmatched[0]["observed_intensity"] == pytest.approx(40.0)
    assert unmatched[0]["predicted_intensity"] is None


def test_ion_mixture_fit_invalid_formula_422(client: TestClient):
    resp = client.post(
        "/api/v1/ion-mixture-fit", json=_payload(formula="C6Xx2")
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["position"] == 2


def test_ion_mixture_fit_peak_validation_422(client: TestClient):
    # not sorted by increasing m/z
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(
            peaks=[
                {"mz": 182.0, "intensity": 1.0},
                {"mz": 181.0, "intensity": 1.0},
            ]
        ),
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_peaks"

    # all-zero intensities
    resp = client.post(
        "/api/v1/ion-mixture-fit",
        json=_payload(peaks=[{"mz": 181.0, "intensity": 0.0}]),
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_peaks"


def test_ion_mixture_fit_parameter_bounds_422(client: TestClient):
    # candidate count must be 2..10
    assert (
        client.post(
            "/api/v1/ion-mixture-fit",
            json=_payload(candidates=[{"adduct": "+H", "charge": 1}]),
        ).status_code
        == 422
    )
    assert (
        client.post(
            "/api/v1/ion-mixture-fit",
            json=_payload(
                candidates=[{"adduct": "+H", "charge": 1}] * 11
            ),
        ).status_code
        == 422
    )

    # step_percent out of 0.5..10
    for step in (0.4, 10.1):
        assert (
            client.post(
                "/api/v1/ion-mixture-fit", json=_payload(step_percent=step)
            ).status_code
            == 422
        )

    # non-positive tolerance
    assert (
        client.post(
            "/api/v1/ion-mixture-fit", json=_payload(mass_tolerance=0.0)
        ).status_code
        == 422
    )

    # zero charge
    assert (
        client.post(
            "/api/v1/ion-mixture-fit",
            json=_payload(
                candidates=[
                    {"adduct": "+H", "charge": 0},
                    {"adduct": "+H", "charge": 1},
                ]
            ),
        ).status_code
        == 422
    )

    # empty peaks / too many peaks
    assert (
        client.post(
            "/api/v1/ion-mixture-fit", json=_payload(peaks=[])
        ).status_code
        == 422
    )
    too_many = [{"mz": float(i), "intensity": 1.0} for i in range(201)]
    assert (
        client.post(
            "/api/v1/ion-mixture-fit", json=_payload(peaks=too_many)
        ).status_code
        == 422
    )
