"""Record filters: QUAL, FILTER column, chromosome region, INFO expression.

Region strings use 1-based inclusive coordinates: 'chr1:100-200'.
INFO expressions: 'KEY<op>value' with op in ==,!=,>=,<=,>,<,= .
Numeric comparison is used when both sides parse as float, else string.
"""

import re

_REGION_RE = re.compile(r"^([^:]+)(?::(\d+)(?:-(\d+))?)?$")
_EXPR_RE = re.compile(r"^([A-Za-z_][\w.]*)\s*(==|!=|>=|<=|>|<|=)\s*(.+)$")


class Region(object):
    __slots__ = ("chrom", "start", "end")

    def __init__(self, chrom, start=None, end=None):
        self.chrom = chrom
        self.start = start
        self.end = end

    @classmethod
    def parse(cls, text):
        m = _REGION_RE.match(text.strip())
        if not m:
            raise ValueError("invalid region %r (want chrom[:start[-end]])" % text)
        chrom, start, end = m.group(1), m.group(2), m.group(3)
        start = int(start) if start else None
        end = int(end) if end else start
        if start is not None and end < start:
            raise ValueError("invalid region %r: end < start" % text)
        return cls(chrom, start, end)

    def contains(self, chrom, pos):
        if chrom != self.chrom:
            return False
        if self.start is None:
            return True
        return self.start <= pos <= self.end


class InfoExpression(object):
    def __init__(self, key, op, raw_value):
        self.key = key
        self.op = op
        self.raw_value = raw_value
        try:
            self.num_value = float(raw_value)
        except ValueError:
            self.num_value = None

    @classmethod
    def parse(cls, text):
        m = _EXPR_RE.match(text.strip())
        if not m:
            raise ValueError("invalid INFO expression %r" % text)
        return cls(m.group(1), m.group(2), m.group(3).strip())

    def matches(self, record):
        value = record.info.get(self.key)
        if value is None:
            return False
        if value is True:  # Flag
            return self.op in ("=", "==") and self.raw_value in ("1", "true", "True")
        if self.num_value is not None:
            try:
                lhs = float(value)
            except (TypeError, ValueError):
                return False
            rhs = self.num_value
        else:
            lhs = str(value)
            rhs = self.raw_value
        if self.op in ("=", "=="):
            return lhs == rhs
        if self.op == "!=":
            return lhs != rhs
        if self.op == ">":
            return lhs > rhs
        if self.op == "<":
            return lhs < rhs
        if self.op == ">=":
            return lhs >= rhs
        if self.op == "<=":
            return lhs <= rhs
        return False


class RecordFilter(object):
    """Composite AND filter over QUAL, FILTER column, regions, INFO exprs."""

    def __init__(self, min_qual=None, require_pass=False, regions=None,
                 info_exprs=None):
        self.min_qual = min_qual
        self.require_pass = require_pass
        self.regions = [Region.parse(r) if isinstance(r, str) else r
                        for r in (regions or [])]
        self.info_exprs = [InfoExpression.parse(e) if isinstance(e, str) else e
                           for e in (info_exprs or [])]

    def accepts(self, record):
        if self.min_qual is not None:
            if record.qual is None or record.qual < self.min_qual:
                return False
        if self.require_pass and not record.is_pass:
            return False
        if self.regions:
            if not any(r.contains(record.chrom, record.pos)
                       for r in self.regions):
                return False
        for expr in self.info_exprs:
            if not expr.matches(record):
                return False
        return True
