import pytest

from vcftool import VCFReader, VCFHeaderError, VCFRecordError
from vcftool.genotype import Genotype


def test_parse_header_and_samples(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        header = reader.header
        assert header.samples == ["HG001", "HG002", "HG003"]
        assert header.contigs["chr1"] == 1000000
        assert header.info_defs["DP"].type == "Integer"
        assert header.info_defs["AF"].number == "A"
        assert header.info_defs["DB"].type == "Flag"
        assert header.format_defs["AD"].number == "R"


def test_parse_snp_record(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        rec = next(reader)
        assert (rec.chrom, rec.pos, rec.ref, rec.alts) == ("chr1", 100, "A", ["G"])
        assert rec.qual == 50.0
        assert rec.filter == "PASS"
        assert rec.info["DP"] == 30
        assert rec.info["AF"] == [0.5]
        assert rec.is_snp and not rec.is_indel and not rec.is_multiallelic
        assert rec.samples[0]["GT"] == Genotype([0, 1])
        assert rec.samples[0]["DP"] == 12
        assert rec.samples[2]["GT"].phased is True


def test_parse_indel_record(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        next(reader)
        rec = next(reader)
        assert rec.ref == "AC" and rec.alts == ["A"]
        assert rec.is_indel and not rec.is_snp
        assert rec.filter == "q10"
        assert rec.samples[1]["GT"].is_missing
        assert rec.samples[1]["DP"] == 0


def test_parse_multiallelic_record(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        recs = list(reader)
        rec = recs[2]
        assert rec.is_multiallelic
        assert rec.alts == ["T", "G"]
        assert rec.info["AF"] == [0.3, 0.1]
        assert rec.info["DB"] is True
        assert rec.samples[1]["GT"] == Genotype([1, 2])
        assert rec.samples[1]["AD"] == [8, 6, 6]


def test_missing_and_phased_genotypes(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        recs = list(reader)
        rec = recs[3]  # chr2:150
        assert rec.qual is None
        assert rec.samples[0]["GT"].phased is True
        assert rec.samples[0]["GT"].is_het
        assert rec.samples[2]["GT"].is_missing


def test_no_sample_vcf(data_path):
    with VCFReader(data_path("no_samples.vcf")) as reader:
        assert reader.header.samples == []
        recs = list(reader)
        assert len(recs) == 2
        assert recs[0].info["DP"] == 5


def test_reader_from_lines(vcf_lines):
    lines = vcf_lines(["chr1\t1\t.\tA\tG\t10\tPASS\tDP=4\tGT\t0/1\t1/1"])
    reader = VCFReader(lines)
    rec = next(iter(reader))
    assert rec.pos == 1
    assert rec.samples[1]["GT"] == Genotype([1, 1])


def test_roundtrip_to_line(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        rec = next(reader)
        line = rec.to_line()
        rec2 = reader.parse_record(line, 999)
        assert rec2.pos == rec.pos
        assert rec2.info["DP"] == 30
        assert str(rec2.samples[0]["GT"]) == "0/1"
