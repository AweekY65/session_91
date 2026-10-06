"""Parsing: SNP, INDEL, multiallelic, phased/unphased, missing genotypes."""

import pytest

from vcftool import Genotype, VCFReader
from vcftool.errors import VCFHeaderError, VCFRecordError


def test_parse_snp_record(vcf_factory):
    rec = "chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30\tGT:DP\t0/1:15\t1/1:15"
    reader = VCFReader(vcf_factory([rec]))
    record = next(reader)
    assert record.chrom == "chr1"
    assert record.pos == 100
    assert record.ref == "A"
    assert record.alts == ["G"]
    assert record.qual == 50.0
    assert record.is_pass
    assert record.variant_type == "SNP"
    assert record.info["DP"] == "30"
    assert record.info_value("DP", int) == 30
    gts = record.genotypes()
    assert gts[0] == Genotype((0, 1), phased=False)
    assert gts[1] == Genotype((1, 1), phased=False)
    assert gts[0].is_heterozygous
    assert not gts[1].is_heterozygous


def test_parse_indel_record(vcf_factory):
    rec = "chr1\t200\t.\tAT\tA\t.\tPASS\tDP=8\tGT\t0/1\t0/0"
    record = next(VCFReader(vcf_factory([rec])))
    assert record.variant_type == "INDEL"
    assert record.qual is None  # '.' QUAL


def test_parse_multiallelic_record(vcf_factory):
    rec = "chr1\t300\t.\tC\tT,G\t99\tPASS\tDP=20;AF=0.25,0.25\tGT\t1/2\t0/2"
    record = next(VCFReader(vcf_factory([rec])))
    assert record.alts == ["T", "G"]
    assert record.n_alleles == 3
    gts = record.genotypes()
    assert gts[0].alleles == (1, 2)
    assert gts[1].alleles == (0, 2)


def test_phased_and_unphased(vcf_factory):
    rec = "chr1\t400\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0|1\t1/0"
    record = next(VCFReader(vcf_factory([rec])))
    gts = record.genotypes()
    assert gts[0].phased and gts[0].alleles == (0, 1)
    assert not gts[1].phased and gts[1].alleles == (1, 0)
    assert str(gts[0]) == "0|1"
    assert str(gts[1]) == "1/0"


def test_missing_genotypes(vcf_factory):
    rec = "chr1\t500\t.\tA\tG\t50\tPASS\tDP=10\tGT\t./.\t.|."
    record = next(VCFReader(vcf_factory([rec])))
    gts = record.genotypes()
    assert gts[0].is_missing and not gts[0].is_called
    assert gts[1].is_missing
    assert gts[1].phased  # separator preserved
    assert gts[0].called_alleles() == []


def test_haploid_and_triploid(vcf_factory):
    rec = "chr1\t600\t.\tA\tG\t50\tPASS\tDP=10\tGT\t1\t0/1/1"
    record = next(VCFReader(vcf_factory([rec])))
    gts = record.genotypes()
    assert gts[0].ploidy == 1 and gts[0].is_called
    assert gts[1].ploidy == 3 and gts[1].is_heterozygous


def test_flag_info_and_empty_filter(vcf_factory):
    rec = "chr1\t700\t.\tA\tG\t50\t.\tDB\tGT\t0/0\t0/0"
    record = next(VCFReader(vcf_factory([rec])))
    assert record.info["DB"] is True
    assert record.filters == []
    assert record.is_pass


def test_header_metadata(vcf_factory):
    reader = VCFReader(vcf_factory([]))
    header = reader.header
    assert header.fileformat == "VCFv4.2"
    assert header.samples == ["S1", "S2"]
    assert header.info["DP"].get("Type") == "Integer"
    assert header.formats["AD"].get("Number") == "R"
    assert "chr1" in header.contigs


def test_missing_column_header_is_error():
    import io
    bad = io.StringIO("##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\t.\t.\n")
    with pytest.raises(VCFHeaderError):
        VCFReader(bad)


def test_streaming_iteration_order(vcf_factory):
    recs = [
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/0",
        "chr1\t2\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/0",
        "chr1\t3\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/0",
    ]
    positions = [r.pos for r in VCFReader(vcf_factory(recs))]
    assert positions == [1, 2, 3]
