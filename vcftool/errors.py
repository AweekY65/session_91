"""Error types for vcftool. All parse errors carry a 1-based line number."""


class VCFError(Exception):
    """Base error for VCF parsing/validation failures."""

    def __init__(self, message, line_number=None):
        self.line_number = line_number
        if line_number is not None:
            message = "line %d: %s" % (line_number, message)
        super().__init__(message)


class VCFHeaderError(VCFError):
    """Malformed or missing VCF header."""


class VCFRecordError(VCFError):
    """Malformed VCF data line (bad columns, POS, FORMAT, truncation...)."""
