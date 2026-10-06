"""Error handling: every malformed input raises with a 1-based line number."""

import pytest

from vcftool import VCFReader
from vcftool.errors import VCFRecordError


def _line_no(vcf_factory, records, **kw):
    reader = VCFReader(vcf_factory(records, **kw))
    with pytest.raises(VCFRecordError) as excinfo:
        for record in reader:
            record.genotypes()  # GT validation happens on genotype access
    return excinfo.value


def test_too_few_columns(vcf_factory):
    err = _line_no(vcf_factory, ["chr1\t100\trs1\tA\tG"])
    assert err.line_no is not None
    assert "columns" in str(err)


def test_truncated_line_missing_sample(vcf_factory):
    # header has 2 samples; record carries only FORMAT + 1 sample
    err = _line_no(
        vcf_factory,
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=5\tGT\t0/1"],
    )
    assert "truncated" in str(err)
    assert err.line_no is not None


def test_truncated_sample_field(vcf_factory):
    # FORMAT declares GT:DP but sample has only one value
    err = _line_no(
        vcf_factory,
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=5\tGT:DP\t0/1\t0/1:9"],
    )
    assert "truncated" in str(err)


def test_illegal_pos_non_integer(vcf_factory):
    err = _line_no(
        vcf_factory,
        ["chr1\tabc\t.\tA\tG\t50\tPASS\tDP=5\tGT\t0/1\t0/1"],
    )
    assert "POS" in str(err)
    assert err.line_no is not None


def test_illegal_pos_zero_and_negative(vcf_factory):
    for pos in ("0", "-3"):
        err = _line_no(
            vcf_factory,
            ["chr1\t%s\t.\tA\tG\t50\tPASS\tDP=5\tGT\t0/1\t0/1" % pos],
        )
        assert "POS" in str(err)


def test_unknown_format_key(vcf_factory):
    err = _line_no(
        vcf_factory,
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=5\tGT:ZZ\t0/1:3\t0/1:3"],
    )
    assert "unknown FORMAT" in str(err)
    assert "'ZZ'" in str(err)


def test_gt_allele_out_of_range(vcf_factory):
    err = _line_no(
        vcf_factory,
        ["chr1\t100\t.\tA\tG\t50\tPASS\tDP=5\tGT\t0/2\t0/1"],
    )
    assert "out of range" in str(err)


def test_line_number_points_at_bad_record(vcf_factory):
    # 11 meta lines (see conftest) + 1 column line = 12; bad record is
    # the 3rd data line -> line 15.
    recs = [
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/0",
        "chr1\t2\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/0",
        "chr1\t3\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0",  # truncated
    ]
    err = _line_no(vcf_factory, recs)
    assert err.line_no == 15
    assert "line 15" in str(err)
