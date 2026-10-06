"""Filtering by QUAL, FILTER column, region and INFO expression."""

import pytest

from vcftool import RecordFilter, VCFReader

RECS = [
    "chr1\t100\t.\tA\tG\t5\tPASS\tDP=30\tGT\t0/1\t0/1",
    "chr1\t200\t.\tA\tG\t50\tq10\tDP=3\tGT\t0/1\t0/1",
    "chr1\t300\t.\tA\tG\t80\tPASS\tDP=25\tGT\t0/1\t0/1",
    "chr2\t100\t.\tA\tG\t90\tPASS\tDP=40\tGT\t0/1\t0/1",
]


def _run(vcf_factory, record_filter):
    reader = VCFReader(vcf_factory(RECS))
    return [(r.chrom, r.pos) for r in reader if record_filter.accepts(r)]


def test_min_qual(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(min_qual=50))
    assert kept == [("chr1", 200), ("chr1", 300), ("chr2", 100)]


def test_require_pass(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(require_pass=True))
    assert kept == [("chr1", 100), ("chr1", 300), ("chr2", 100)]


def test_region(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(regions=["chr1:150-350"]))
    assert kept == [("chr1", 200), ("chr1", 300)]


def test_region_single_position_and_whole_chrom(vcf_factory):
    assert _run(vcf_factory, RecordFilter(regions=["chr1:200"])) == [("chr1", 200)]
    assert _run(vcf_factory, RecordFilter(regions=["chr2"])) == [("chr2", 100)]


def test_multiple_regions(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(regions=["chr1:100", "chr2:1-150"]))
    assert kept == [("chr1", 100), ("chr2", 100)]


def test_info_expression_numeric(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(info_exprs=["DP>=25"]))
    assert kept == [("chr1", 100), ("chr1", 300), ("chr2", 100)]


def test_info_expression_string(vcf_factory):
    recs = ["chr1\t1\t.\tA\tG\t10\tPASS\tTYPE=SNP;DP=5\tGT\t0/1\t0/1",
            "chr1\t2\t.\tAT\tA\t10\tPASS\tTYPE=DEL;DP=5\tGT\t0/1\t0/1"]
    reader = VCFReader(vcf_factory(recs))
    kept = [r.pos for r in reader
            if RecordFilter(info_exprs=["TYPE=SNP"]).accepts(r)]
    assert kept == [1]


def test_combined_filters(vcf_factory):
    kept = _run(vcf_factory, RecordFilter(
        min_qual=50, require_pass=True, regions=["chr1:1-1000"],
        info_exprs=["DP>=20"]))
    assert kept == [("chr1", 300)]


def test_invalid_region_and_expression():
    with pytest.raises(ValueError):
        RecordFilter(regions=["chr1:50-10"])
    with pytest.raises(ValueError):
        RecordFilter(info_exprs=["not an expression"])
