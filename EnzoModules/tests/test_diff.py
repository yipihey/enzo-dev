import math

from enzomodules.diff import BITWISE, Tolerance, compare, isclose


def test_bitwise():
    assert isclose(1.0, 1.0, BITWISE)
    assert not isclose(1.0, math.nextafter(1.0, 2.0), BITWISE)


def test_relative():
    assert isclose(1.0, 1.0 + 1e-13, Tolerance(rtol=1e-12))
    assert not isclose(1.0, 1.1, Tolerance(rtol=1e-12))


def test_nan_roundtrips_equal():
    assert isclose(math.nan, math.nan, BITWISE)
    assert not isclose(math.nan, 1.0, BITWISE)


def test_compare_pass_and_fail():
    a = [1.0, 2.0, 3.0]
    assert compare(a, a, BITWISE).ok
    r = compare([1.0, 2.0, 3.0], [1.0, 2.0, 3.5], Tolerance(rtol=1e-9))
    assert not r.ok
    assert r.nfail == 1
    assert r.worst == 3
    assert math.isclose(r.maxabs, 0.5)


def test_compareresult_is_truthy():
    assert compare([1.0], [1.0], BITWISE)        # __bool__
    assert not compare([1.0], [2.0], BITWISE)
