"""Streaming VCF parser. Pure local: reads only file objects or local paths.

Validates, with 1-based line numbers in every error:
  * wrong column counts (too few / truncated / too many)
  * illegal POS (non-integer or < 1)
  * FORMAT keys not declared in the header (GT is always allowed)
  * truncated sample fields (fewer values than FORMAT keys)
  * GT allele indices out of range for the record's alleles
"""

from .errors import VCFHeaderError, VCFRecordError
from .genotype import Genotype

FIXED_COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


class MetaDeclaration(object):
    """A structured ##KEY=<ID=...,Number=...,Type=...,Description=...> line."""

    __slots__ = ("key", "fields")

    def __init__(self, key, fields):
        self.key = key
        self.fields = fields

    def get(self, name, default=None):
        return self.fields.get(name, default)


def _parse_structured_meta(line, line_no):
    """Parse '##INFO=<ID=DP,Number=1,Type=Integer,Description="...">'."""
    key, _, rest = line[2:].partition("=")
    rest = rest.strip()
    if not (rest.startswith("<") and rest.endswith(">")):
        return MetaDeclaration(key, {})
    body = rest[1:-1]
    fields = {}
    token, current, in_quotes = [], [], False
    parts = []
    for ch in body:
        if ch == '"':
            in_quotes = not in_quotes
            current.append(ch)
        elif ch == "," and not in_quotes:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    parts.append("".join(current))
    for part in parts:
        k, _, v = part.partition("=")
        fields[k.strip()] = v.strip().strip('"')
    return MetaDeclaration(key, fields)


class VCFHeader(object):
    def __init__(self):
        self.meta_lines = []          # raw '##' lines, order preserved
        self.fileformat = None
        self.info = {}                # ID -> MetaDeclaration
        self.formats = {}             # ID -> MetaDeclaration
        self.filters = {}             # ID -> MetaDeclaration
        self.contigs = {}             # ID -> MetaDeclaration
        self.samples = []             # sample names from the #CHROM line
        self.column_line = None       # raw '#CHROM ...' line

    def add_meta_line(self, line, line_no):
        self.meta_lines.append(line)
        if line.startswith("##fileformat="):
            self.fileformat = line.split("=", 1)[1].strip()
            return
        for prefix, table in (
            ("##INFO=<", self.info),
            ("##FORMAT=<", self.formats),
            ("##FILTER=<", self.filters),
            ("##contig=<", self.contigs),
        ):
            if line.startswith(prefix):
                decl = _parse_structured_meta(line, line_no)
                decl_id = decl.get("ID")
                if decl_id:
                    table[decl_id] = decl
                return

    def set_column_line(self, line, line_no):
        self.column_line = line
        cols = line.lstrip("#").rstrip("\n").split("\t")
        if len(cols) < 8 or tuple(cols[:8]) != FIXED_COLUMNS:
            raise VCFHeaderError(
                "column header must start with #CHROM POS ID REF ALT QUAL "
                "FILTER INFO, got %r" % line.rstrip("\n"),
                line_no,
            )
        if len(cols) > 8:
            if cols[8] != "FORMAT":
                raise VCFHeaderError(
                    "9th column must be FORMAT when samples are present",
                    line_no,
                )
            self.samples = cols[9:]

    @property
    def n_samples(self):
        return len(self.samples)


class VCFRecord(object):
    __slots__ = (
        "chrom", "pos", "vid", "ref", "alts", "qual", "filters", "info",
        "format_keys", "sample_data", "line_no",
    )

    def __init__(self, chrom, pos, vid, ref, alts, qual, filters, info,
                 format_keys, sample_data, line_no):
        self.chrom = chrom
        self.pos = pos
        self.vid = vid
        self.ref = ref
        self.alts = alts
        self.qual = qual            # float or None when '.'
        self.filters = filters      # list of str; [] when '.'
        self.info = info            # dict: key -> str value or True (Flag)
        self.format_keys = format_keys
        self.sample_data = sample_data  # list of dicts, one per sample
        self.line_no = line_no

    @property
    def n_alleles(self):
        return 1 + len(self.alts)

    @property
    def is_pass(self):
        return not self.filters or self.filters == ["PASS"]

    @property
    def variant_type(self):
        if any(a.startswith("<") for a in self.alts):
            return "SYMBOLIC"
        if all(len(self.ref) == 1 and len(a) == 1 for a in self.alts):
            return "SNP"
        return "INDEL"

    def genotypes(self):
        """Parsed GT per sample, in header sample order. None if no FORMAT."""
        out = []
        for data in self.sample_data:
            gt_text = data.get("GT")
            if gt_text is None:
                out.append(None)
            else:
                try:
                    gt = Genotype.parse(gt_text)
                except ValueError as exc:
                    raise VCFRecordError(str(exc), self.line_no)
                for allele in gt.called_alleles():
                    if allele >= self.n_alleles:
                        raise VCFRecordError(
                            "GT allele index %d out of range: record has "
                            "%d allele(s) (REF + %d ALT)"
                            % (allele, self.n_alleles, len(self.alts)),
                            self.line_no,
                        )
                out.append(gt)
        return out

    def info_value(self, key, cast=None):
        """Typed INFO accessor. Returns None when absent."""
        value = self.info.get(key)
        if value is None or value is True:
            return value
        if cast is not None:
            try:
                return cast(value)
            except (TypeError, ValueError):
                return None
        return value


class VCFReader(object):
    """Iterate VCFRecord objects from a path or text file object.

    Streaming: records are parsed lazily, one line at a time, so files of
    any size can be processed with constant memory.
    """

    def __init__(self, source):
        self._owns_file = False
        if isinstance(source, str):
            self._fh = open(source, "r")
            self._owns_file = True
        else:
            self._fh = source
        self.header = VCFHeader()
        self._line_no = 0
        self._read_header()

    def _read_header(self):
        for raw in self._fh:
            self._line_no += 1
            line = raw.rstrip("\n")
            if line.startswith("##"):
                self.header.add_meta_line(line, self._line_no)
            elif line.startswith("#"):
                self.header.set_column_line(line, self._line_no)
                return
            elif line.strip() == "":
                continue
            else:
                raise VCFHeaderError(
                    "expected header line starting with '#', got %r" % line,
                    self._line_no,
                )
        raise VCFHeaderError("missing #CHROM column header line")

    def __iter__(self):
        return self

    def __next__(self):
        for raw in self._fh:
            self._line_no += 1
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            return self._parse_record(line, self._line_no)
        if self._owns_file:
            self._fh.close()
        raise StopIteration

    next = __next__

    def _parse_record(self, line, line_no):
        cols = line.split("\t")
        if len(cols) < 8:
            raise VCFRecordError(
                "expected at least 8 tab-separated columns, got %d "
                "(truncated or malformed line)" % len(cols),
                line_no,
            )
        expected = 8 + (1 + self.header.n_samples if self.header.n_samples else 0)
        if self.header.n_samples:
            if len(cols) < expected:
                raise VCFRecordError(
                    "truncated line: expected %d columns (8 fixed + FORMAT + "
                    "%d samples), got %d"
                    % (expected, self.header.n_samples, len(cols)),
                    line_no,
                )
            if len(cols) > expected:
                raise VCFRecordError(
                    "too many columns: expected %d, got %d"
                    % (expected, len(cols)),
                    line_no,
                )
        elif len(cols) > 8:
            raise VCFRecordError(
                "header declares no samples but data line has %d columns"
                % len(cols),
                line_no,
            )

        chrom = cols[0]
        try:
            pos = int(cols[1])
        except ValueError:
            raise VCFRecordError(
                "illegal POS %r: must be a positive integer" % cols[1],
                line_no,
            )
        if pos < 1:
            raise VCFRecordError(
                "illegal POS %d: must be >= 1" % pos, line_no,
            )

        vid = cols[2]
        ref = cols[3]
        if not ref or ref == ".":
            raise VCFRecordError("REF allele is missing", line_no)
        alts = [] if cols[4] in (".", "") else cols[4].split(",")

        qual = None if cols[5] == "." else None
        if cols[5] != ".":
            try:
                qual = float(cols[5])
            except ValueError:
                raise VCFRecordError(
                    "illegal QUAL %r: must be a number or '.'" % cols[5],
                    line_no,
                )

        filters = [] if cols[6] in (".", "") else cols[6].split(";")

        info = {}
        if cols[7] not in (".", ""):
            for item in cols[7].split(";"):
                key, eq, value = item.partition("=")
                if not key:
                    raise VCFRecordError("empty INFO key", line_no)
                info[key] = value if eq else True

        format_keys = []
        sample_data = []
        if self.header.n_samples:
            format_keys = cols[8].split(":") if cols[8] else []
            for key in format_keys:
                if key != "GT" and key not in self.header.formats:
                    raise VCFRecordError(
                        "unknown FORMAT key %r: not declared in any "
                        "##FORMAT header line" % key,
                        line_no,
                    )
            for sample_idx in range(self.header.n_samples):
                raw_values = cols[9 + sample_idx].split(":")
                if len(raw_values) < len(format_keys):
                    raise VCFRecordError(
                        "truncated sample field for sample %r: FORMAT "
                        "declares %d keys but only %d value(s) present"
                        % (self.header.samples[sample_idx],
                           len(format_keys), len(raw_values)),
                        line_no,
                    )
                if len(raw_values) > len(format_keys):
                    raise VCFRecordError(
                        "sample %r has %d value(s) but FORMAT declares "
                        "only %d key(s)"
                        % (self.header.samples[sample_idx],
                           len(raw_values), len(format_keys)),
                        line_no,
                    )
                sample_data.append(dict(zip(format_keys, raw_values)))

        return VCFRecord(
            chrom, pos, vid, ref, alts, qual, filters, info,
            format_keys, sample_data, line_no,
        )

    def close(self):
        if self._owns_file and not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
