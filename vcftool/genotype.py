"""Genotype (GT) semantics.

A genotype is a list of allele indices (ints) or None for a missing allele,
plus a `phased` flag: True when alleles were separated by '|', False for '/'.
Mixed separators in one GT are treated as unphased.
A GT is "missing" when every allele is None (e.g. '.', './.', '.|.').
"""

MISSING = "."


class Genotype(object):
    __slots__ = ("alleles", "phased")

    def __init__(self, alleles, phased=False):
        self.alleles = tuple(alleles)
        self.phased = phased

    @classmethod
    def parse(cls, text, line_number=None):
        from .errors import VCFRecordError

        if text in ("", MISSING):
            return cls([None], phased=False)
        phased = "|" in text and "/" not in text
        parts = text.replace("|", "/").split("/")
        alleles = []
        for part in parts:
            if part == MISSING:
                alleles.append(None)
            else:
                try:
                    allele = int(part)
                except ValueError:
                    raise VCFRecordError(
                        "invalid GT allele %r" % part, line_number)
                if allele < 0:
                    raise VCFRecordError(
                        "negative GT allele %d" % allele, line_number)
                alleles.append(allele)
        return cls(alleles, phased=phased)

    @property
    def ploidy(self):
        return len(self.alleles)

    @property
    def is_missing(self):
        return all(a is None for a in self.alleles)

    @property
    def is_called(self):
        return not self.is_missing

    @property
    def is_het(self):
        called = [a for a in self.alleles if a is not None]
        return len(set(called)) > 1

    def __str__(self):
        sep = "|" if self.phased else "/"
        return sep.join(MISSING if a is None else str(a) for a in self.alleles)

    def __eq__(self, other):
        return (isinstance(other, Genotype)
                and self.alleles == other.alleles
                and self.phased == other.phased)

    def __repr__(self):
        return "Genotype(%s)" % str(self)
