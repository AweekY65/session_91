import pytest

from vcftool.errors import VCFRecordError
from vcftool.genotype import Genotype


def test_unphased():
    gt = Genotype.parse("0/1")
    assert gt.alleles == (0, 1)
    assert not gt.phased
    assert gt.is_het and gt.is_called


def test_phased():
    gt = Genotype.parse("1|0")
    assert gt.alleles == (1, 0)
    assert gt.phased
    assert gt.is_het


def test_missing_forms():
    for text in (".", "./.", ".|."):
        gt = Genotype.parse(text)
        assert gt.is_missing
        assert not gt.is_called


def test_partially_missing():
    gt = Genotype.parse("0/.")
    assert gt.is_called
    assert not gt.is_het
    assert gt.alleles == (0, None)


def test_multiallelic():
    gt = Genotype.parse("1/2")
    assert gt.is_het
    assert gt.alleles == (1, 2)


def test_haploid():
    gt = Genotype.parse("1")
    assert gt.ploidy == 1
    assert not gt.is_het


def test_str_roundtrip():
    assert str(Genotype.parse("0|1")) == "0|1"
    assert str(Genotype.parse("./.")) == "./."


def test_invalid_allele():
    with pytest.raises(VCFRecordError):
        Genotype.parse("0/x", line_number=7)


def test_error_has_line_number():
    with pytest.raises(VCFRecordError) as exc:
        Genotype.parse("0/x", line_number=7)
    assert "line 7" in str(exc.value)
