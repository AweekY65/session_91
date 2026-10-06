import pytest

from vcftool import VCFReader, VCFHeaderError, VCFRecordError


def read_all(lines):
    return list(VCFReader(lines))


def test_wrong_column_count(vcf_lines):
    lines = vcf_lines(["chr1\t1\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "wrong column count" in str(exc.value)
    assert "line 8" in str(exc.value)


def test_truncated_line(vcf_lines):
    lines = vcf_lines(["chr1\t1\t.\tA"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "truncated line" in str(exc.value)
    assert "line 8" in str(exc.value)


def test_invalid_pos_not_integer(vcf_lines):
    lines = vcf_lines(["chr1\tabc\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "invalid POS" in str(exc.value)
    assert "line 8" in str(exc.value)


def test_invalid_pos_zero(vcf_lines):
    lines = vcf_lines(["chr1\t0\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "invalid POS 0" in str(exc.value)


def test_unknown_format_key(vcf_lines):
    lines = vcf_lines(["chr1\t1\t.\tA\tG\t10\tPASS\tDP=4\tGT:ZZ\t0/1:3\t0/0:4"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "unknown FORMAT key 'ZZ'" in str(exc.value)
    assert "line 8" in str(exc.value)


def test_invalid_ref(vcf_lines):
    lines = vcf_lines(["chr1\t5\t.\tAX\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0"])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "invalid REF" in str(exc.value)


def test_missing_chrom_header():
    with pytest.raises(VCFHeaderError):
        read_all(["##fileformat=VCFv4.2\n"])


def test_malformed_meta_line():
    with pytest.raises(VCFHeaderError) as exc:
        read_all(["##badmetalinex\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"])
    assert "line 1" in str(exc.value)


def test_line_numbers_track_multiple_records(vcf_lines):
    lines = vcf_lines([
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0",
        "chr1\t2\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0",
        "chr1\tbad\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t0/0",
    ])
    with pytest.raises(VCFRecordError) as exc:
        read_all(lines)
    assert "line 10" in str(exc.value)
