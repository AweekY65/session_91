"""Error types for vcftool. All parse errors carry a 1-based line number."""


class VCFError(Exception):
    """Base error for all vcftool failures."""

    def __init__(self, message, line_no=None):
        self.line_no = line_no
        if line_no is not None:
            message = "line %d: %s" % (line_no, message)
        super().__init__(message)


class VCFHeaderError(VCFError):
    """Malformed or missing VCF header."""


class VCFRecordError(VCFError):
    """Malformed VCF data line (bad columns, POS, FORMAT, truncation)."""
