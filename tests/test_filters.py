from vcftool import VCFReader, RecordFilter, parse_region, eval_info_expression
from vcftool.errors import VCFError
import pytest


def load(path):
    with VCFReader(path) as reader:
        return list(reader)


def test_min_qual(data_path):
    records = load(data_path("basic.vcf"))
    kept = list(RecordFilter(min_qual=40).filter(records))
    assert [r.pos for r in kept] == [100, 300]


def test_min_qual_excludes_missing_qual(data_path):
    records = load(data_path("basic.vcf"))
    kept = list(RecordFilter(min_qual=0).filter(records))
    assert all(r.qual is not None for r in kept)


def test_keep_filter(data_path):
    records = load(data_path("basic.vcf"))
    kept = list(RecordFilter(keep_filters=["PASS"]).filter(records))
    assert [r.pos for r in kept] == [100, 300, 150]
    kept = list(RecordFilter(keep_filters=["q10"]).filter(records))
    assert [r.pos for r in kept] == [200]


def test_region(data_path):
    records = load(data_path("basic.vcf"))
    kept = list(RecordFilter(regions=["chr1:150-400"]).filter(records))
    assert [r.pos for r in kept] == [200, 300]
    kept = list(RecordFilter(regions=["chr2"]).filter(records))
    assert [r.pos for r in kept] == [150, 250]


def test_region_overlaps_indel_span(data_path):
    records = load(data_path("basic.vcf"))
    # indel at chr1:200 REF=AC spans 200-201; region starting at 201 hits it
    kept = list(RecordFilter(regions=["chr1:201-201"]).filter(records))
    assert [r.pos for r in kept] == [200]


def test_info_expressions(data_path):
    records = load(data_path("basic.vcf"))
    kept = list(RecordFilter(info_expressions=["DP>=30"]).filter(records))
    assert [r.pos for r in kept] == [100, 300]
    kept = list(RecordFilter(info_expressions=["TYPE=SNP"]).filter(records))
    assert [r.pos for r in kept] == [100, 300, 250]
    kept = list(RecordFilter(info_expressions=["DB"]).filter(records))
    assert [r.pos for r in kept] == [300]
    kept = list(RecordFilter(info_expressions=["AF<0.2"]).filter(records))
    assert [r.pos for r in kept] == [300]  # AF=0.3,0.1 -> any element < 0.2


def test_combined_filters(data_path):
    records = load(data_path("basic.vcf"))
    filt = RecordFilter(min_qual=10, regions=["chr1"], info_expressions=["DP>20"])
    kept = list(filt.filter(records))
    assert [r.pos for r in kept] == [100, 200, 300]


def test_parse_region_errors():
    with pytest.raises(VCFError):
        parse_region("chr1:0-10")
    with pytest.raises(VCFError):
        parse_region("chr1:20-10")
    with pytest.raises(VCFError):
        parse_region("chr1:abc")
    region = parse_region("chr1:1,000-2,000")
    assert (region.start, region.end) == (1000, 2000)


def test_eval_info_expression_ne(data_path):
    records = load(data_path("basic.vcf"))
    rec = records[0]
    assert eval_info_expression(rec, "TYPE!=INDEL")
    assert not eval_info_expression(rec, "MISSING_KEY>1")
