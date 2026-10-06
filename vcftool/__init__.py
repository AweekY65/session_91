"""vcftool: local-only VCF parsing, filtering, splitting and statistics."""

from .errors import VCFError, VCFHeaderError, VCFRecordError
from .genotype import Genotype
from .parser import VCFReader, VCFHeader, VariantRecord, open_vcf
from .filters import RecordFilter, parse_region, eval_info_expression
from .split import split_record, split_multiallelic
from .stats import StatsCollector, collect_stats, record_allele_counts

__version__ = "0.1.0"

__all__ = [
    "VCFError", "VCFHeaderError", "VCFRecordError",
    "Genotype", "VCFReader", "VCFHeader", "VariantRecord", "open_vcf",
    "RecordFilter", "parse_region", "eval_info_expression",
    "split_record", "split_multiallelic",
    "StatsCollector", "collect_stats", "record_allele_counts",
]
