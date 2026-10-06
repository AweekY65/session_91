"""Streaming behaviour on large / unbounded inputs."""

import itertools

from vcftool import StatsCollector, VCFReader

from conftest import HEADER


class InfiniteVCF(object):
    """File-like object yielding an endless valid VCF, one line at a time."""

    def __init__(self):
        header = HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        records = (
            "chr1\t%d\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n" % (i + 1)
            for i in itertools.count()
        )
        self._lines = itertools.chain(header.splitlines(keepends=True), records)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._lines)


def test_reader_is_lazy_on_unbounded_input():
    reader = VCFReader(InfiniteVCF())
    stream = iter(reader)
    first_ten = [next(stream).pos for _ in range(10)]
    assert first_ten == list(range(1, 11))


def test_large_generated_vcf(tmp_path):
    n_records = 50000
    path = tmp_path / "large.vcf"
    with open(path, "w") as fh:
        fh.write(HEADER)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
        for i in range(n_records):
            gt1 = "0/1" if i % 2 == 0 else "1/1"
            gt2 = "./." if i % 5 == 0 else "0/0"
            fh.write("chr1\t%d\t.\tA\tG\t60\tPASS\tDP=20\tGT\t%s\t%s\n"
                     % (i + 1, gt1, gt2))

    collector = None
    with VCFReader(str(path)) as reader:
        collector = StatsCollector(reader.header.samples)
        for record in reader:
            collector.add_record(record)

    assert collector.n_records == n_records
    s1, s2 = collector.samples["S1"], collector.samples["S2"]
    assert s1.called == n_records
    assert s1.call_rate == 1.0
    # S2 missing every 5th record
    assert s2.called == n_records - n_records // 5
    assert s2.call_rate == 0.8
    # S1 het on even i (0/1): half of all records
    assert s1.het == n_records // 2
    # AN over all sites: S1 gives 2 alleles/record, S2 2 alleles unless missing
    total_an = sum(site.an for site in collector.site_stats)
    assert total_an == 2 * n_records + 2 * (n_records - n_records // 5)
