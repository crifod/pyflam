"""Tests for the OpenFOAM dictionary serializer (no OpenFOAM)."""

from __future__ import annotations

import pytest

from pyflam.cfd import dictwriter as dw


def test_scalars_and_bools():
    out = dw.render({"a": 1, "b": 2.5, "flag": True, "off": False})
    assert "a               1;" in out
    assert "b               2.5;" in out
    assert "flag            true;" in out
    assert "off             false;" in out


def test_list_and_verbatim_string():
    out = dw.render({"g": [0, 0, -9.81], "val": "uniform (1 2 3)"})
    assert "g               (0 0 -9.81);" in out
    assert "val             uniform (1 2 3);" in out  # string emitted verbatim


def test_nested_dict_block():
    out = dw.render({"SIMPLE": {"pRefCell": 0, "pRefValue": 0}})
    assert "SIMPLE\n{" in out
    assert "    pRefCell        0;" in out


def test_foam_file_header_and_footer():
    text = dw.foam_file("volScalarField", "k", {"dimensions": "[0 2 -2 0 0 0 0]"},
                        location="0")
    assert "FoamFile" in text
    assert "class       volScalarField;" in text
    assert 'location    "0";' in text
    assert "object      k;" in text
    assert text.rstrip().endswith("//")


def test_unserializable_raises():
    with pytest.raises(TypeError):
        dw.render({"x": object()})
