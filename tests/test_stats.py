"""Statistics: call rate, heterozygosity, allele count, allele frequency."""

import pytest

from vcftool import StatsCollector, VCFReader


def _collect(vcf_factory, records, samples=("S1", "S2")):
    reader = VCFReader(vcf_factory(records, samples=samples))
    collector = StatsCollector(reader.header.samples)
    for record in reader:
        collector.add_record(record)
    return collector


def test_call_rate_with_missing(vcf_factory):
    recs = [
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/1\t./.",
        "chr1\t2\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/0\t0/1",
        "chr1\t3\t.\tA\tG\t10\tPASS\tDP=1\tGT\t./.\t0/1",
        "chr1\t4\t.\tA\tG\t10\tPASS\tDP=1\tGT\t1/1\t0/0",
    ]
    c = _collect(vcf_factory, recs)
    s1, s2 = c.samples["S1"], c.samples["S2"]
    assert s1.total == 4 and s1.called == 3
    assert s1.call_rate == pytest.approx(0.75)
    assert s2.called == 3 and s2.call_rate == pytest.approx(0.75)


def test_heterozygosity(vcf_factory):
    recs = [
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/1\t0/0",
        "chr1\t2\t.\tA\tG\t10\tPASS\tDP=1\tGT\t1/1\t0|1",
        "chr1\t3\t.\tA\tG\t10\tPASS\tDP=1\tGT\t./.\t0/0",
    ]
    c = _collect(vcf_factory, recs)
    s1, s2 = c.samples["S1"], c.samples["S2"]
    # S1: 2 called, 1 het -> 0.5 ; S2: 3 called, 1 het (phased counts)
    assert s1.heterozygosity == pytest.approx(0.5)
    assert s2.heterozygosity == pytest.approx(1 / 3)


def test_allele_count_and_frequency(vcf_factory):
    recs = [
        "chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/1\t1/1",
        "chr1\t2\t.\tC\tT\t10\tPASS\tDP=1\tGT\t0/0\t./.",
    ]
    c = _collect(vcf_factory, recs)
    site1, site2 = c.site_stats
    # site1: alleles 0,1,1,1 -> AC=[1,3], AN=4, AF=[0.25,0.75]
    assert site1.ac == [1, 3]
    assert site1.an == 4
    assert site1.allele_frequency == pytest.approx([0.25, 0.75])
    # site2: only 0/0 called -> AC=[2,0], AN=2, AF=[1.0,0.0]
    assert site2.ac == [2, 0]
    assert site2.an == 2
    assert site2.allele_frequency == pytest.approx([1.0, 0.0])


def test_multiallelic_allele_counts(vcf_factory):
    recs = ["chr1\t1\t.\tA\tG,T\t10\tPASS\tDP=1\tGT\t1/2\t0/2"]
    c = _collect(vcf_factory, recs)
    site = c.site_stats[0]
    assert site.ac == [1, 1, 2]
    assert site.an == 4
    assert site.allele_frequency == pytest.approx([0.25, 0.25, 0.5])


def test_haploid_contributes_one_allele(vcf_factory):
    recs = ["chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t1\t0/0"]
    c = _collect(vcf_factory, recs)
    site = c.site_stats[0]
    assert site.ac == [2, 1]
    assert site.an == 3


def test_summary_output_contains_stats(vcf_factory):
    recs = ["chr1\t1\t.\tA\tG\t10\tPASS\tDP=1\tGT\t0/1\t1/1"]
    c = _collect(vcf_factory, recs)
    text = "\n".join(c.summary_lines())
    assert "call_rate" in text
    assert "S1\t1.0000" in text
    assert "chr1:1" in text
