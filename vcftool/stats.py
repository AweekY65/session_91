"""Variant statistics.

Definitions
-----------
Per-sample (over all records in the stream):
  call_rate      = called_genotypes / total_records
  heterozygosity = heterozygous_calls / called_genotypes   (diploid het only)

Per-record (allele index 0 = REF, 1..n = ALT):
  AC (allele count)      = number of non-missing GT alleles equal to index
  AN (allele number)     = total number of non-missing GT alleles
  AF (allele frequency)  = AC / AN  (None when AN == 0)
"""


class SampleStats(object):
    __slots__ = ("name", "called", "het", "total")

    def __init__(self, name):
        self.name = name
        self.called = 0
        self.het = 0
        self.total = 0

    @property
    def call_rate(self):
        return self.called / self.total if self.total else 0.0

    @property
    def heterozygosity(self):
        return self.het / self.called if self.called else 0.0


class RecordStats(object):
    __slots__ = ("record", "ac", "an")

    def __init__(self, record, ac, an):
        self.record = record
        self.ac = ac  # list of counts, index 0 = REF
        self.an = an

    @property
    def af(self):
        if self.an == 0:
            return [None] * len(self.ac)
        return [count / self.an for count in self.ac]


def record_allele_counts(record):
    """Return (ac, an) for one record."""
    n_alleles = 1 + len(record.alts)
    ac = [0] * n_alleles
    an = 0
    for sample in record.samples:
        gt = sample.get("GT")
        if gt is None:
            continue
        for allele in gt.alleles:
            if allele is None:
                continue
            if allele < n_alleles:
                ac[allele] += 1
            an += 1
    return ac, an


class StatsCollector(object):
    """Streaming accumulator over a record iterator."""

    def __init__(self, sample_names):
        self.sample_stats = [SampleStats(n) for n in sample_names]
        self.n_records = 0
        self.n_snps = 0
        self.n_indels = 0
        self.n_multiallelic = 0
        self.record_stats = []  # per-record RecordStats, in stream order

    def add(self, record):
        self.n_records += 1
        if record.is_snp:
            self.n_snps += 1
        elif record.is_indel:
            self.n_indels += 1
        if record.is_multiallelic:
            self.n_multiallelic += 1
        for stats, sample in zip(self.sample_stats, record.samples):
            stats.total += 1
            gt = sample.get("GT")
            if gt is not None and gt.is_called:
                stats.called += 1
                if gt.ploidy == 2 and gt.is_het:
                    stats.het += 1
        ac, an = record_allele_counts(record)
        self.record_stats.append(RecordStats(record, ac, an))
        return self.record_stats[-1]


def collect_stats(records, sample_names):
    collector = StatsCollector(sample_names)
    for record in records:
        collector.add(record)
    return collector


def format_report(collector):
    """Render a plain-text report of collected statistics."""
    lines = []
    lines.append("== Variant summary ==")
    lines.append("records:        %d" % collector.n_records)
    lines.append("SNPs:           %d" % collector.n_snps)
    lines.append("INDELs:         %d" % collector.n_indels)
    lines.append("multiallelic:   %d" % collector.n_multiallelic)
    lines.append("")
    if collector.sample_stats:
        lines.append("== Per-sample stats ==")
        lines.append("%-16s %8s %8s %8s %10s" % (
            "sample", "called", "total", "callrate", "het"))
        for s in collector.sample_stats:
            lines.append("%-16s %8d %8d %8.4f %10.4f" % (
                s.name, s.called, s.total, s.call_rate, s.heterozygosity))
        lines.append("")
    lines.append("== Per-record allele stats ==")
    header = "%-10s %10s  %-5s %-5s  %s" % ("CHROM", "POS", "REF", "ALT",
                                            "AC(0..n)/AN -> AF")
    lines.append(header)
    for rs in collector.record_stats:
        rec = rs.record
        acs = ",".join(str(c) for c in rs.ac)
        afs = ",".join("NA" if f is None else "%.4f" % f for f in rs.af)
        lines.append("%-10s %10d  %-5s %-5s  %s/%d -> %s" % (
            rec.chrom, rec.pos, rec.ref, ",".join(rec.alts) or ".",
            acs, rs.an, afs))
    return "\n".join(lines)
