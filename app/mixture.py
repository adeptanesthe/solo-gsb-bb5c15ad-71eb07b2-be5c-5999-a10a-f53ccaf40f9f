"""Two-component ion-mixture fitting over candidate pairs and mixing ratios.

Given the same neutral formula, each valid candidate (adduct, charge) is
turned into a base-peak-normalised m/z envelope with the same machinery as
``/api/v1/ion-identify`` (:mod:`app.identify`). Every unordered pair of
distinct valid candidates is then mixed over a 0-100% grid for the first
candidate's share: both envelopes keep their own base-peak normalisation
(base peak = 100) while intensities are scaled by the mixing fractions,
within-tolerance peaks of the superposition are fused, and the mixture is
base-peak normalised again. The mixed spectrum is paired one-to-one with the
observed (likewise base-peak-normalised) peaks using the existing
:func:`app.enrichment.pair_and_score`, so unmatched peaks on either side
count at zero intensity in the MSE.

For each candidate pair the grid point with the smallest MSE wins (a tie
within floating-point noise keeps the lower first-candidate fraction, as the
grid is visited ascending). Pairs are ranked by MSE with the lexicographically
smaller candidate index pair as tie-break. The near-equality tolerance only
absorbs float drift from the intensity-weighted peak fusion (superposing two
identical envelopes at different fractions cannot change the spectrum in
exact arithmetic); real MSE differences are orders of magnitude larger.
"""

from __future__ import annotations

from dataclasses import dataclass

from .enrichment import abundance_grid, normalise_base_peak, pair_and_score
from .envelope import merge_peaks
from .identify import IonCandidatePrediction

# (m/z, intensity) peak on the m/z axis, mirroring identify.MzPeak.
MzPeak = tuple[float, float]


@dataclass
class MixtureScore:
    """Best grid-point score for one unordered pair of candidates."""

    first_index: int
    second_index: int
    # Share (percent) of the first candidate in the winning mixture.
    first_fraction_percent: float
    mse: float
    # (observed_index | None, predicted_index | None) pairs at the best point
    pairs: list[tuple[int | None, int | None]]
    # Base-peak renormalised mixed m/z spectrum at the best grid point
    predicted: list[MzPeak]


# Relative MSE tolerance below which two grid points are regarded as tied
# (superposing identical envelopes can jitter the intensity-weighted merge
# means by floating-point noise at different mixing fractions).
MSE_TIE_REL_TOLERANCE = 1e-12


def _is_strict_better(mse: float, best: float) -> bool:
    """True when ``mse`` beats ``best`` beyond floating-point tie noise.

    A difference smaller than 1e-12 of the larger magnitude is an exact tie
    for ranking purposes; ties keep the grid point / candidate pair visited
    first (lower fraction, then smaller index pair).
    """
    return mse < best - MSE_TIE_REL_TOLERANCE * max(1.0, abs(best), abs(mse))


def mix_envelopes(
    peaks_a: list[MzPeak],
    peaks_b: list[MzPeak],
    fraction_a_percent: float,
    tolerance: float,
) -> list[MzPeak]:
    """Superpose two base-peak-normalised m/z envelopes.

    Each envelope keeps its own base-peak (= 100) normalisation; intensities
    are scaled by the candidate's fraction (``fraction_a_percent/100`` for
    ``peaks_a`` and the complement for ``peaks_b``). The scaled peaks are
    fused when closer than ``tolerance`` (m/z) with the existing
    :func:`merge_peaks` rule (intensity-weighted mean mass, summed
    intensity), and the superposition is base-peak normalised again.
    """
    fraction_a = fraction_a_percent / 100.0
    fraction_b = 1.0 - fraction_a
    combined: list[MzPeak] = [
        (mz, intensity * fraction_a)
        for mz, intensity in peaks_a
        if fraction_a > 0.0
    ]
    if fraction_b > 0.0:
        combined.extend(
            (mz, intensity * fraction_b) for mz, intensity in peaks_b
        )
    # At grid endpoints one fraction is exactly 0; that side contributes no
    # peaks, so the single-component mixture is simply the other envelope.
    if not combined:
        return []
    merged = merge_peaks(combined, tolerance)
    return normalise_base_peak(merged)


def fit_ion_mixture(
    predictions: list[tuple[int, IonCandidatePrediction]],
    observed: list[MzPeak],
    step_percent: float,
    mass_tolerance: float,
) -> tuple[list[MixtureScore], int]:
    """Score every unordered pair of valid candidates over the mixing grid.

    ``predictions`` are ``(request index, prediction)`` pairs in ascending
    request-index order; ``observed`` must already be base-peak normalised.
    Returns the per-pair best scores ranked by (MSE, first index, second
    index) plus the number of grid points scanned for each pair.
    """
    grid = abundance_grid(step_percent)

    results: list[MixtureScore] = []
    # Ascending (first_index, second_index) visit order realises the index
    # tie-break; within a pair the ascending grid realises the fraction one.
    for k, (index_a, prediction_a) in enumerate(predictions):
        for index_b, prediction_b in predictions[k + 1 :]:
            best: MixtureScore | None = None
            for fraction_percent in grid:
                mixed = mix_envelopes(
                    prediction_a.peaks,
                    prediction_b.peaks,
                    fraction_percent,
                    mass_tolerance,
                )
                pairs, mse = pair_and_score(observed, mixed, mass_tolerance)
                if best is None or _is_strict_better(mse, best.mse):
                    best = MixtureScore(
                        index_a,
                        index_b,
                        fraction_percent,
                        mse,
                        pairs,
                        mixed,
                    )
            # Every pair visits at least the 0% grid point, so best exists.
            assert best is not None
            results.append(best)

    # Rank by MSE, treating differences below the relative tie tolerance as
    # equal; within a tie block the visit order (ascending index pair) is
    # kept. Implemented by numbering ascending-MSE tie groups, then stably
    # re-sorting the visit-ordered results by group number.
    by_mse = sorted(results, key=lambda r: r.mse)
    group_of: dict[tuple[int, int], int] = {}
    group = -1
    group_min_mse = 0.0
    for score in by_mse:
        if group < 0 or _is_strict_better(score.mse, group_min_mse):
            group += 1
            group_min_mse = score.mse
        group_of[(score.first_index, score.second_index)] = group

    ranked = sorted(
        results,
        key=lambda r: group_of[(r.first_index, r.second_index)],
    )
    return ranked, len(grid)


__all__ = [
    "MixtureScore",
    "fit_ion_mixture",
    "mix_envelopes",
]
