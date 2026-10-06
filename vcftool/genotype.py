"""Genotype model: phased/unphased, missing, multiallelic.

A genotype is a tuple of allele indices (0 = REF, 1..n = ALT) or None for
a missing allele, plus a `phased` flag (True for '0|1', False for '0/1').
"""

_GT_SEP_PHASED = "|"
_GT_SEP_UNPHASED = "/"


class Genotype(object):
    __slots__ = ("alleles", "phased")

    def __init__(self, alleles, phased=False):
        self.alleles = tuple(alleles)
        self.phased = phased

    @classmethod
    def parse(cls, text):
        """Parse a GT string such as '0/1', '1|0', './.', '.', '2'."""
        if text in (".", ""):
            return cls((None,), False)
        if _GT_SEP_PHASED in text:
            phased = True
            parts = text.split(_GT_SEP_PHASED)
        else:
            phased = False
            parts = text.split(_GT_SEP_UNPHASED)
        alleles = []
        for part in parts:
            if part == ".":
                alleles.append(None)
            else:
                try:
                    idx = int(part)
                except ValueError:
                    raise ValueError("invalid GT allele %r" % part)
                if idx < 0:
                    raise ValueError("invalid GT allele %r" % part)
                alleles.append(idx)
        return cls(alleles, phased)

    @property
    def ploidy(self):
        return len(self.alleles)

    @property
    def is_called(self):
        """A genotype is called when it has no missing allele."""
        return bool(self.alleles) and all(a is not None for a in self.alleles)

    @property
    def is_missing(self):
        return not self.is_called

    @property
    def is_heterozygous(self):
        """Called genotype with at least two distinct allele indices."""
        if not self.is_called:
            return False
        return len(set(self.alleles)) > 1

    def called_alleles(self):
        return [a for a in self.alleles if a is not None]

    def __str__(self):
        sep = _GT_SEP_PHASED if self.phased else _GT_SEP_UNPHASED
        return sep.join("." if a is None else str(a) for a in self.alleles)

    def __eq__(self, other):
        return (
            isinstance(other, Genotype)
            and self.alleles == other.alleles
            and self.phased == other.phased
        )

    def __repr__(self):
        return "Genotype(%s)" % str(self)
