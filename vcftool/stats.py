"""Per-sample and per-site statistics computed over a stream of records.

Definitions
-----------
call rate(sample)   = called genotypes / total records
heterozygosity(s)   = heterozygous called genotypes / called genotypes
AC(site, allele)    = count of that allele over all called alleles
AN(site)            = total number of called (non-missing) alleles
AF(site, allele)    = AC / AN
"""

from collections import OrderedDict


class SampleStats(object):
    __slots__ = ("name", "total", "called", "het")

    def __init__(self, name):
        self.name = name
        self.total = 0
        self.called = 0
        self.het = 0

    @property
    def call_rate(self):
        return self.called / self.total if self.total else 0.0

    @property
    def heterozygosity(self):
        return self.het / self.called if self.called else 0.0


class SiteStats(object):
    __slots__ = ("chrom", "pos", "vid", "ref", "alts", "an", "ac")

    def __init__(self, record):
        self.chrom = record.chrom
        self.pos = record.pos
        self.vid = record.vid
        self.ref = record.ref
        self.alts = list(record.alts)
        self.an = 0
        self.ac = [0] * record.n_alleles

    def observe(self, genotype):
        for allele in genotype.called_alleles():
            self.ac[allele] += 1
            self.an += 1

    @property
    def allele_frequency(self):
        if self.an == 0:
            return [0.0] * len(self.ac)
        return [count / self.an for count in self.ac]


class StatsCollector(object):
    """Accumulates sample and site statistics over records, one at a time."""

    def __init__(self, sample_names):
        self.samples = OrderedDict((n, SampleStats(n)) for n in sample_names)
        self.n_records = 0
        self.site_stats = []

    def add_record(self, record):
        genotypes = record.genotypes()
        self.n_records += 1
        site = SiteStats(record)
        for idx, sample in enumerate(self.samples.values()):
            sample.total += 1
            gt = genotypes[idx] if idx < len(genotypes) else None
            if gt is None or gt.is_missing:
                continue
            sample.called += 1
            if gt.is_heterozygous:
                sample.het += 1
            site.observe(gt)
        self.site_stats.append(site)

    def summary_lines(self):
        lines = []
        lines.append("records\t%d" % self.n_records)
        lines.append("")
        lines.append("sample\tcall_rate\tcalled\ttotal\thet\theterozygosity")
        for sample in self.samples.values():
            lines.append("%s\t%.4f\t%d\t%d\t%d\t%.4f" % (
                sample.name, sample.call_rate, sample.called, sample.total,
                sample.het, sample.heterozygosity,
            ))
        lines.append("")
        lines.append("site\tAN\tAC\tAF")
        for site in self.site_stats:
            name = "%s:%d" % (site.chrom, site.pos)
            ac = ",".join(str(c) for c in site.ac)
            af = ",".join("%.4f" % f for f in site.allele_frequency)
            lines.append("%s\t%d\t%s\t%s" % (name, site.an, ac, af))
        return lines
