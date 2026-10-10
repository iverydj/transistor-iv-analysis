import numpy as np
import pytest

from transistor_iv import CurveError, monotonic_segments, select_segment
from transistor_iv.sweep import on_direction


def test_segments_share_turning_points():
    v = [0, -2, -4, -2, 0, 2, 0]
    assert [s.tolist() for s in monotonic_segments(v)] == [[0, 1, 2], [2, 3, 4, 5], [5, 6]]


def test_repeated_point_separates_segments():
    v = [0, 1, 2, 2, 1, 0]
    assert [s.tolist() for s in monotonic_segments(v)] == [[0, 1, 2], [3, 4, 5]]


def test_select_longest_first_of_equal_and_out_of_range():
    v = np.r_[0, -1, np.arange(-2, 3), np.arange(1, -3, -1)]   # 0 -> -2 -> 2 -> -2
    segs = monotonic_segments(v)
    assert len(segs) == 3 and len(segs[1]) == len(segs[2])
    np.testing.assert_array_equal(select_segment(v, -1), segs[1])
    with pytest.raises(CurveError, match="segment 3"):
        select_segment(v, 3)


def test_on_direction_mirrors_p_type():
    v = np.array([-2.0, -1.0, 0.0])
    i = np.array([-1e-6, -1e-9, -1e-12])
    u, a, s = on_direction(v, i, "p")
    assert s == -1
    np.testing.assert_array_equal(u, [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(a, [1e-12, 1e-9, 1e-6])
    with pytest.raises(CurveError, match="monotonic"):
        on_direction([0, 1, 0], [1, 2, 3], "n")
