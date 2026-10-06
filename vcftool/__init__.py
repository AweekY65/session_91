"""vcftool: local-only VCF parsing, filtering, splitting and statistics.

No network access of any kind: all input comes from local files or
in-memory file objects.
"""

from .errors import VCFError, VCFHeaderError, VCFRecordError
from .genotype import Genotype
from .parser import VCFReader, VCFRecord, VCFHeader
from .stats import StatsCollector
from .filters import RecordFilter, Region, InfoExpression
from .split import split_multiallelic, format_record

__version__ = "0.1.0"

__all__ = [
    "VCFError", "VCFHeaderError", "VCFRecordError",
    "Genotype", "VCFReader", "VCFRecord", "VCFHeader",
    "StatsCollector", "RecordFilter", "Region", "InfoExpression",
    "split_multiallelic", "format_record",
]
