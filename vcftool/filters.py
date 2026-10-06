"""Record filters: QUAL, FILTER column, chromosome region, INFO expressions.

Region strings use 1-based inclusive coordinates: "chr1:100-200" or "chr1".
INFO expressions: "KEY" (presence), "KEY=value", "KEY!=value",
"KEY>v", "KEY>=v", "KEY<v", "KEY<=v" (numeric comparison when possible).
"""

import re

from .errors import VCFError

_REGION_RE = re.compile(r"^(?P<chrom>[^:]+)(:(?P<start>[\d,]+)(-(?P<end>[\d,]+))?)?$")
_EXPR_RE = re.compile(r"^(?P<key>[^!<>=]+?)(?P<op>>=|<=|!=|>|<|=)(?P<value>.*)$")


class Region(object):
    def __init__(self, chrom, start=None, end=None):
        self.chrom = chrom
        self.start = start
        self.end = end

    def contains(self, record):
        if record.chrom != self.chrom:
            return False
        if self.start is not None and record.end < self.start:
            return False
        if self.end is not None and record.pos > self.end:
            return False
        return True


def parse_region(text):
    match = _REGION_RE.match(text.strip())
    if not match:
        raise VCFError("invalid region %r (expected chr[:start[-end]])" % text)
    start = match.group("start")
    end = match.group("end")
    start = int(start.replace(",", "")) if start else None
    end = int(end.replace(",", "")) if end else None
    if start is not None and start < 1:
        raise VCFError("invalid region %r: start must be >= 1" % text)
    if end is not None and start is not None and end < start:
        raise VCFError("invalid region %r: end < start" % text)
    return Region(match.group("chrom"), start, end)


def _coerce_pair(value, target):
    """Try to compare numerically; fall back to string comparison."""
    try:
        return float(value), float(target)
    except (TypeError, ValueError):
        return str(value), str(target)


def eval_info_expression(record, expr):
    expr = expr.strip()
    match = _EXPR_RE.match(expr)
    if not match:
        # bare key: presence test
        return expr in record.info
    key, op, raw = match.group("key"), match.group("op"), match.group("value")
    if key not in record.info:
        return False
    value = record.info[key]
    if isinstance(value, list):
        return any(_compare(v, op, raw) for v in value if v is not None)
    return _compare(value, op, raw)


def _compare(value, op, raw):
    left, right = _coerce_pair(value, raw)
    if op == "=":
        return left == right
    if op == "!=":
        return left != right
    if op == ">":
        return left > right
    if op == ">=":
        return left >= right
    if op == "<":
        return left < right
    if op == "<=":
        return left <= right
    raise VCFError("unsupported operator %r" % op)  # pragma: no cover


class RecordFilter(object):
    """Composable predicate over VariantRecord."""

    def __init__(self, min_qual=None, keep_filters=None, regions=None,
                 info_expressions=None):
        self.min_qual = min_qual
        # keep_filters: iterable of FILTER values to accept; "PASS" and "."
        # both mean "passing".
        self.keep_filters = (set(keep_filters) if keep_filters is not None
                             else None)
        self.regions = [parse_region(r) for r in (regions or [])]
        self.info_expressions = list(info_expressions or [])

    def accept(self, record):
        if self.min_qual is not None:
            if record.qual is None or record.qual < self.min_qual:
                return False
        if self.keep_filters is not None:
            filt = record.filter
            passing = filt in ("PASS", ".")
            if not (filt in self.keep_filters
                    or (passing and "PASS" in self.keep_filters)):
                return False
        if self.regions:
            if not any(r.contains(record) for r in self.regions):
                return False
        for expr in self.info_expressions:
            if not eval_info_expression(record, expr):
                return False
        return True

    def filter(self, records):
        for record in records:
            if self.accept(record):
                yield record
