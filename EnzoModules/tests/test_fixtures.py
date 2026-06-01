import os

import enzomodules
from enzomodules.fixtures import load_dir, load_fixture, save_fixture


def _fixture_dir():
    return enzomodules.fixturedir("Hydro", "twoshock")


def test_fixtures_present_and_parsed():
    fixes = load_dir(_fixture_dir())
    assert fixes, "no fixtures found"
    by_name = {f.name: f for f in fixes}
    assert "sod" in by_name
    sod = by_name["sod"]
    assert sod["gamma"] == 1.4
    assert sod["dls"] == [1.0]
    assert sod["drs"] == [0.125]
    assert "pbar" in sod and "ubar" in sod


def test_scalar_typing():
    sod = {f.name: f for f in load_dir(_fixture_dir())}["sod"]
    assert isinstance(sod["idim"], int)
    assert isinstance(sod["gamma"], float)
    assert isinstance(sod.scalars["name"], str)


def test_save_load_roundtrip(tmp_path):
    sod = {f.name: f for f in load_dir(_fixture_dir())}["sod"]
    p = os.path.join(tmp_path, "sod.fixture")
    save_fixture(p, sod)
    again = load_fixture(p)
    assert again["pbar"] == sod["pbar"]
    assert again["ubar"] == sod["ubar"]
    assert again["gamma"] == sod["gamma"]
    assert again.name == "sod"
