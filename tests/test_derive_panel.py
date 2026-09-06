"""Peak-calling logic for `pyeuk derive-panel`.

These exercise the pure-numpy core (peak calling, adaptive threshold, BAM sampling)
without pysam or real BAMs, so they run in the base install.
"""

import numpy as np

from pyeuk.amplicon.derive_panel import _sample_bams, call_peaks


def _cov(length, peaks):
    """Build a coverage array of `length` with (start, end, depth) peaks on a 1x background."""
    a = np.ones(length, dtype=np.int64)
    for s, e, d in peaks:
        a[s:e] = d
    return a


def test_calls_two_clear_amplicons():
    cov = {"c1": _cov(1000, [(100, 300, 5000), (600, 800, 8000)])}
    lengths = {"c1": 1000}
    peaks, thresh, gmax = call_peaks(cov, min_depth_frac=0.01, min_abs_depth=30,
                                     min_length=80, merge_gap=25, pad=0, lengths=lengths)
    assert gmax == 8000
    assert thresh == 80  # 0.01 * 8000
    assert len(peaks) == 2
    # sorted by depth desc: the 8000x peak first
    assert (peaks[0]["contig"], peaks[0]["start"], peaks[0]["end"]) == ("c1", 600, 800)
    assert peaks[0]["max_depth"] == 8000
    assert peaks[1]["start"] == 100 and peaks[1]["end"] == 300


def test_threshold_is_data_adaptive():
    # same shape, 10x deeper -> threshold scales, peak count unchanged (no retuning).
    # both depths kept above the absolute floor so the fraction, not the floor, sets it.
    shallow = {"c1": _cov(1000, [(100, 300, 5000)])}
    deep = {"c1": _cov(1000, [(100, 300, 50000)])}
    lengths = {"c1": 1000}
    p1, t1, _ = call_peaks(shallow, 0.01, 30, 80, 25, 0, lengths)
    p2, t2, _ = call_peaks(deep, 0.01, 30, 80, 25, 0, lengths)
    assert len(p1) == len(p2) == 1
    assert t1 == 50 and t2 == 500  # threshold tracks depth (10x), no retuning


def test_absolute_floor_wins_when_fraction_is_tiny():
    cov = {"c1": _cov(500, [(100, 300, 200)])}
    lengths = {"c1": 500}
    _, thresh, _ = call_peaks(cov, min_depth_frac=0.001, min_abs_depth=30, min_length=80,
                              merge_gap=25, pad=0, lengths=lengths)
    assert thresh == 30  # 0.001*200 = 0.2 -> floored at 30


def test_merge_gap_joins_adjacent_amplicons():
    # two peaks 10 bp apart merge into one (an amplicon that dips mid-peak)
    cov = {"c1": _cov(1000, [(100, 250, 5000), (260, 400, 5000)])}
    lengths = {"c1": 1000}
    peaks, _, _ = call_peaks(cov, 0.01, 30, 80, merge_gap=25, pad=0, lengths=lengths)
    assert len(peaks) == 1
    assert peaks[0]["start"] == 100 and peaks[0]["end"] == 400


def test_short_specks_dropped_by_min_length():
    cov = {"c1": _cov(1000, [(100, 300, 5000), (700, 730, 9000)])}  # 2nd is only 30 bp
    lengths = {"c1": 1000}
    peaks, _, _ = call_peaks(cov, 0.01, 30, min_length=80, merge_gap=25, pad=0, lengths=lengths)
    assert len(peaks) == 1
    assert peaks[0]["start"] == 100


def test_pad_extends_and_clamps_to_contig():
    cov = {"c1": _cov(1000, [(5, 200, 5000)])}
    lengths = {"c1": 1000}
    peaks, _, _ = call_peaks(cov, 0.01, 30, 80, 25, pad=50, lengths=lengths)
    assert peaks[0]["start"] == 0        # 5-50 clamped to 0
    assert peaks[0]["end"] == 250        # 200+50


def test_sample_bams_even_stride():
    bams = [f"b{i}" for i in range(100)]
    picked = _sample_bams(bams, 10)
    assert len(picked) == 10
    assert picked[0] == "b0" and picked[1] == "b10"  # strided, not the first 10
    assert _sample_bams(bams, 0) == bams             # 0 = all
    assert _sample_bams(["x", "y"], 10) == ["x", "y"]  # fewer than k -> all
