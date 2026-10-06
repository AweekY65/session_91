import math

from vcftool import VCFReader, collect_stats
from vcftool.stats import format_report, record_allele_counts


def collect(path):
    with VCFReader(path) as reader:
        return collect_stats(iter(reader), reader.header.samples)


def test_sample_call_rate(data_path):
    stats = collect(data_path("basic.vcf"))
    by_name = {s.name: s for s in stats.sample_stats}
    # HG002 has a missing GT at chr1:200 -> 4/5 called
    assert math.isclose(by_name["HG002"].call_rate, 4 / 5)
    assert math.isclose(by_name["HG001"].call_rate, 1.0)
    # HG003 has a missing GT at chr2:150 -> 4/5 called
    assert math.isclose(by_name["HG003"].call_rate, 4 / 5)


def test_sample_heterozygosity(data_path):
    stats = collect(data_path("basic.vcf"))
    by_name = {s.name: s for s in stats.sample_stats}
    # HG001: calls 0/1, 0/1, 0/1, 0|1, 0/0 -> 4 het of 5 called
    assert math.isclose(by_name["HG001"].heterozygosity, 4 / 5)
    # HG003: 0|0, 0/0, 2/2, ./., 0/1 -> 1 het of 4 called
    assert math.isclose(by_name["HG003"].heterozygosity, 1 / 4)


def test_allele_counts_and_frequency(data_path):
    stats = collect(data_path("basic.vcf"))
    rs = stats.record_stats[0]  # chr1:100 A>G, GTs 0/1, 1/1, 0|0
    assert rs.ac == [3, 3]
    assert rs.an == 6
    assert rs.af == [0.5, 0.5]


def test_allele_counts_multiallelic(data_path):
    stats = collect(data_path("basic.vcf"))
    rs = stats.record_stats[2]  # chr1:300 C>T,G: 0/1, 1/2, 2/2
    assert rs.ac == [1, 2, 3]
    assert rs.an == 6
    assert math.isclose(rs.af[2], 0.5)


def test_allele_counts_missing_excluded(data_path):
    stats = collect(data_path("basic.vcf"))
    rs = stats.record_stats[3]  # chr2:150 G>GA: 0|1, 1|1, ./.
    assert rs.ac == [1, 3]
    assert rs.an == 4
    assert math.isclose(rs.af[1], 0.75)


def test_variant_type_counts(data_path):
    stats = collect(data_path("basic.vcf"))
    assert stats.n_records == 5
    assert stats.n_snps == 3
    assert stats.n_indels == 2
    assert stats.n_multiallelic == 1


def test_report_renders(data_path):
    stats = collect(data_path("basic.vcf"))
    report = format_report(stats)
    assert "== Variant summary ==" in report
    assert "HG001" in report
    assert "callrate" in report


def test_record_allele_counts_direct(data_path):
    with VCFReader(data_path("basic.vcf")) as reader:
        rec = next(reader)
        ac, an = record_allele_counts(rec)
        assert ac == [3, 3] and an == 6
