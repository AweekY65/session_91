"""Streaming behaviour over a generated, larger VCF (no GUI, no network)."""

import math
import tracemalloc

from vcftool import VCFReader, RecordFilter, collect_stats
from vcftool.split import split_multiallelic

N_RECORDS = 20000


def generate_vcf_lines(n=N_RECORDS):
    """Deterministically generate a large in-memory VCF."""
    yield "##fileformat=VCFv4.2\n"
    yield "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
    yield "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    yield "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    gts = ["0/0", "0/1", "1/1", "./."]
    for i in range(n):
        if i % 10 == 9:
            ref, alt = "A", "G,T"  # multiallelic every 10th record
        elif i % 7 == 6:
            ref, alt = "AC", "A"   # indel
        else:
            ref, alt = "A", "G"
        qual = (i % 100) + 1
        filt = "PASS" if i % 5 else "q10"
        gt1 = gts[i % 4]
        gt2 = gts[(i + 1) % 4]
        yield "chr1\t%d\t.\t%s\t%s\t%d\t%s\tDP=%d\tGT\t%s\t%s\n" % (
            i + 1, ref, alt, qual, filt, i % 50, gt1, gt2)


def test_streaming_parse_and_stats():
    reader = VCFReader(generate_vcf_lines())
    collector = collect_stats(iter(reader), reader.header.samples)
    assert collector.n_records == N_RECORDS
    assert collector.n_multiallelic == N_RECORDS // 10
    s1, s2 = collector.sample_stats
    # GT pattern cycles with period 4; './.' is gts[3]
    expected_called_s1 = sum(1 for i in range(N_RECORDS) if i % 4 != 3)
    expected_called_s2 = sum(1 for i in range(N_RECORDS) if (i + 1) % 4 != 3)
    assert s1.called == expected_called_s1
    assert s2.called == expected_called_s2
    assert math.isclose(s1.call_rate, expected_called_s1 / N_RECORDS)


def test_streaming_memory_bounded():
    """Peak memory must not grow with record count: process records lazily."""
    tracemalloc.start()
    reader = VCFReader(generate_vcf_lines())
    filt = RecordFilter(min_qual=50)
    n_kept = 0
    an_sum = 0
    for record in filt.filter(reader):
        n_kept += 1
        for sample in record.samples:
            gt = sample.get("GT")
            if gt and gt.is_called:
                an_sum += gt.ploidy
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert n_kept > 0
    # 20k records would take many MB if buffered; lazy processing stays small
    assert peak < 5 * 1024 * 1024


def test_streaming_split():
    reader = VCFReader(generate_vcf_lines())
    n_out = sum(1 for _ in split_multiallelic(iter(reader)))
    expected = N_RECORDS + (N_RECORDS // 10)  # each triallelic adds one
    assert n_out == expected
