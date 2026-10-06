"""Streaming VCF parser (VCFv4.2 subset).

Everything is read from local file objects or iterables of strings; no
network access of any kind is performed.
"""

import re

from .errors import VCFHeaderError, VCFRecordError
from .genotype import Genotype

_META_RE = re.compile(r"^##(?P<key>[^=]+)=(?P<value>.*)$")
_STRUCT_RE = re.compile(r"^<(?P<body>.*)>$")

VALID_BASES = set("ACGTNacgtn")
MIN_COLUMNS = 8  # CHROM POS ID REF ALT QUAL FILTER INFO


class FieldDef(object):
    """Definition of an INFO/FORMAT field from a ##INFO / ##FORMAT header."""

    __slots__ = ("id", "number", "type", "description")

    def __init__(self, id, number, type, description=""):
        self.id = id
        self.number = number  # int, or one of 'A', 'R', 'G', '.'
        self.type = type      # Integer, Float, Flag, Character, String
        self.description = description

    def __repr__(self):
        return "FieldDef(%s, Number=%s, Type=%s)" % (
            self.id, self.number, self.type)


def _parse_structured(value):
    """Parse '<ID=xx,Number=xx,Type=xx,Description="...">' into a dict."""
    match = _STRUCT_RE.match(value.strip())
    if not match:
        return None
    body = match.group("body")
    parts, current, in_quote = [], [], False
    for ch in body:
        if ch == '"':
            in_quote = not in_quote
        elif ch == "," and not in_quote:
            parts.append("".join(current))
            current = []
            continue
        current.append(ch)
    parts.append("".join(current))
    result = {}
    for part in parts:
        if "=" in part:
            key, val = part.split("=", 1)
            result[key.strip()] = val.strip().strip('"')
    return result


class VCFHeader(object):
    def __init__(self):
        self.meta = {}            # key -> list of raw values
        self.info_defs = {}       # ID -> FieldDef
        self.format_defs = {}     # ID -> FieldDef
        self.contigs = {}         # ID -> length or None
        self.samples = []
        self.column_line = None   # raw #CHROM line
        self.line_number = None   # line of the #CHROM header

    def add_meta(self, key, value, line_number):
        self.meta.setdefault(key, []).append(value)
        if key in ("INFO", "FORMAT", "contig"):
            fields = _parse_structured(value)
            if fields is None or "ID" not in fields:
                raise VCFHeaderError(
                    "malformed ##%s declaration: %r" % (key, value),
                    line_number)
            if key == "INFO":
                self.info_defs[fields["ID"]] = FieldDef(
                    fields["ID"], fields.get("Number", "."),
                    fields.get("Type", "String"),
                    fields.get("Description", ""))
            elif key == "FORMAT":
                self.format_defs[fields["ID"]] = FieldDef(
                    fields["ID"], fields.get("Number", "."),
                    fields.get("Type", "String"),
                    fields.get("Description", ""))
            else:
                length = fields.get("length")
                self.contigs[fields["ID"]] = (
                    int(length) if length and length.isdigit() else None)

    @property
    def expected_columns(self):
        return MIN_COLUMNS + (1 + len(self.samples) if self.samples else 0)

    def format_def(self, key):
        # GT is always implicitly known (VCF spec requires it declared, but
        # be lenient for the one universally-defined key).
        if key == "GT":
            return self.format_defs.get(
                key, FieldDef("GT", "1", "String", "Genotype"))
        return self.format_defs.get(key)


def _convert_scalar(text, ftype, line_number, field_id):
    if text == ".":
        return None
    if ftype == "Integer":
        try:
            return int(text)
        except ValueError:
            raise VCFRecordError(
                "INFO/FORMAT field %s: expected Integer, got %r"
                % (field_id, text), line_number)
    if ftype == "Float":
        try:
            return float(text)
        except ValueError:
            raise VCFRecordError(
                "INFO/FORMAT field %s: expected Float, got %r"
                % (field_id, text), line_number)
    return text


def _convert_value(text, fdef, line_number):
    """Convert a raw INFO/FORMAT value using its header definition."""
    if fdef is None:
        return text
    if fdef.type == "Flag":
        return True
    if fdef.number in ("A", "R", "G") or "," in text:
        return [_convert_scalar(p, fdef.type, line_number, fdef.id)
                for p in text.split(",")]
    return _convert_scalar(text, fdef.type, line_number, fdef.id)


class VariantRecord(object):
    """One parsed VCF data line."""

    __slots__ = ("chrom", "pos", "id", "ref", "alts", "qual", "filter",
                 "info", "format_keys", "samples", "line_number", "header")

    def __init__(self, chrom, pos, id, ref, alts, qual, filter, info,
                 format_keys, samples, line_number, header):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alts = alts
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format_keys = format_keys
        self.samples = samples  # list of dicts: FORMAT key -> value
        self.line_number = line_number
        self.header = header

    @property
    def alleles(self):
        return [self.ref] + list(self.alts)

    @property
    def is_snp(self):
        return (len(self.ref) == 1
                and all(len(a) == 1 and a.upper() in VALID_BASES
                        for a in self.alts))

    @property
    def is_indel(self):
        return (len(self.ref) != 1
                or any(len(a) != 1 for a in self.alts)) and not self.is_snp

    @property
    def is_multiallelic(self):
        return len(self.alts) > 1

    @property
    def end(self):
        return self.pos + len(self.ref) - 1

    def genotype(self, sample_name):
        idx = self.header.samples.index(sample_name)
        return self.samples[idx].get("GT")

    def to_line(self):
        cols = [self.chrom, str(self.pos), self.id or ".", self.ref,
                ",".join(self.alts) if self.alts else ".",
                "." if self.qual is None else ("%g" % self.qual),
                self.filter or ".",
                _format_info(self.info)]
        if self.format_keys:
            cols.append(":".join(self.format_keys))
            for sample in self.samples:
                cols.append(":".join(
                    _format_value(sample.get(k)) for k in self.format_keys))
        return "\t".join(cols)


def _format_value(value):
    if value is None:
        return "."
    if isinstance(value, bool):
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join("." if v is None else str(v) for v in value)
    if isinstance(value, Genotype):
        return str(value)
    return str(value)


def _format_info(info):
    if not info:
        return "."
    parts = []
    for key, value in info.items():
        if value is True:
            parts.append(key)
        else:
            parts.append("%s=%s" % (key, _format_value(value)))
    return ";".join(parts)


def parse_info(text, header, line_number):
    info = {}
    if text in ("", "."):
        return info
    for entry in text.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, raw = entry.split("=", 1)
            fdef = header.info_defs.get(key)
            info[key] = _convert_value(raw, fdef, line_number)
        else:
            fdef = header.info_defs.get(entry)
            if fdef is not None and fdef.type == "Flag":
                info[entry] = True
            else:
                info[entry] = True  # tolerate undeclared flags as present
    return info


class VCFReader(object):
    """Streaming reader. Yields VariantRecord one at a time; memory use is
    O(1) in the number of records."""

    def __init__(self, source):
        """`source` is a path (str) or an iterable of text lines."""
        if isinstance(source, str):
            self._stream = open(source, "r")
            self._owns_stream = True
        else:
            self._stream = iter(source)
            self._owns_stream = False
        self.header = VCFHeader()
        self._line_number = 0
        self._pending = None
        self._read_header()

    def _next_line(self):
        for raw in self._stream:
            self._line_number += 1
            return raw.rstrip("\n").rstrip("\r")
        return None

    def _read_header(self):
        saw_column_line = False
        while True:
            line = self._next_line()
            if line is None:
                raise VCFHeaderError(
                    "unexpected end of input: missing #CHROM header line")
            if line.startswith("##"):
                match = _META_RE.match(line)
                if not match:
                    raise VCFHeaderError("malformed meta line: %r" % line,
                                         self._line_number)
                self.header.add_meta(match.group("key"), match.group("value"),
                                     self._line_number)
            elif line.startswith("#CHROM"):
                columns = line.lstrip("#").split("\t")
                if len(columns) < MIN_COLUMNS:
                    raise VCFHeaderError(
                        "column header has %d fields, expected at least %d"
                        % (len(columns), MIN_COLUMNS), self._line_number)
                if len(columns) > MIN_COLUMNS:
                    if columns[MIN_COLUMNS] != "FORMAT":
                        raise VCFHeaderError(
                            "9th column must be FORMAT, got %r"
                            % columns[MIN_COLUMNS], self._line_number)
                    self.header.samples = columns[MIN_COLUMNS + 1:]
                self.header.column_line = line
                self.header.line_number = self._line_number
                saw_column_line = True
                break
            else:
                raise VCFHeaderError(
                    "expected ##meta or #CHROM line, got %r" % line[:60],
                    self._line_number)
        if not saw_column_line:  # pragma: no cover - loop always breaks/raises
            raise VCFHeaderError("missing #CHROM header line")

    def __iter__(self):
        return self

    def __next__(self):
        line = self._next_line()
        if line is None:
            self.close()
            raise StopIteration
        return self.parse_record(line, self._line_number)

    next = __next__

    def close(self):
        if self._owns_stream and not self._stream.closed:
            self._stream.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def parse_record(self, line, line_number):
        if not line.strip():
            raise VCFRecordError("empty line", line_number)
        cols = line.split("\t")
        expected = self.header.expected_columns
        if len(cols) < MIN_COLUMNS:
            raise VCFRecordError(
                "truncated line: %d columns, expected at least %d"
                % (len(cols), MIN_COLUMNS), line_number)
        if self.header.samples and len(cols) != expected:
            raise VCFRecordError(
                "wrong column count: got %d, expected %d (8 fixed + FORMAT "
                "+ %d samples)" % (len(cols), expected,
                                   len(self.header.samples)), line_number)

        chrom = cols[0]
        if not chrom:
            raise VCFRecordError("empty CHROM", line_number)
        try:
            pos = int(cols[1])
        except ValueError:
            raise VCFRecordError("invalid POS %r: not an integer" % cols[1],
                                 line_number)
        if pos < 1:
            raise VCFRecordError("invalid POS %d: must be >= 1" % pos,
                                 line_number)

        ref = cols[3]
        if not ref or ref == ".":
            raise VCFRecordError("missing REF allele", line_number)
        if any(b not in VALID_BASES for b in ref):
            raise VCFRecordError("invalid REF bases %r" % ref, line_number)

        alts = [] if cols[4] in ("", ".") else cols[4].split(",")

        qual = None
        if cols[5] not in ("", "."):
            try:
                qual = float(cols[5])
            except ValueError:
                raise VCFRecordError("invalid QUAL %r" % cols[5], line_number)

        filt = cols[6] if cols[6] not in ("",) else "."
        info = parse_info(cols[7], self.header, line_number)

        format_keys = []
        samples = []
        if len(cols) > MIN_COLUMNS:
            format_keys = cols[MIN_COLUMNS].split(":")
            for key in format_keys:
                if self.header.format_def(key) is None:
                    raise VCFRecordError(
                        "unknown FORMAT key %r: not declared in header" % key,
                        line_number)
            sample_cols = cols[MIN_COLUMNS + 1:]
            if len(sample_cols) != len(self.header.samples):
                raise VCFRecordError(
                    "wrong column count: got %d sample columns, expected %d"
                    % (len(sample_cols), len(self.header.samples)),
                    line_number)
            for raw in sample_cols:
                values = raw.split(":")
                if len(values) > len(format_keys):
                    raise VCFRecordError(
                        "sample column %r has more values than FORMAT keys"
                        % raw, line_number)
                sample = {}
                for key, value in zip(format_keys, values):
                    if key == "GT":
                        sample[key] = Genotype.parse(value, line_number)
                    else:
                        fdef = self.header.format_def(key)
                        sample[key] = _convert_value(value, fdef,
                                                     line_number)
                samples.append(sample)

        return VariantRecord(chrom, pos, cols[2], ref, alts, qual, filt,
                             info, format_keys, samples, line_number,
                             self.header)


def open_vcf(source):
    return VCFReader(source)
