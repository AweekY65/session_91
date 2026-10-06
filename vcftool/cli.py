"""Command line interface: python -m vcftool <stats|filter|split> ..."""

import argparse
import sys

from .errors import VCFError
from .filters import RecordFilter
from .parser import VCFReader
from .split import split_multiallelic
from .stats import StatsCollector, format_report


def _build_filter(args):
    return RecordFilter(
        min_qual=args.min_qual,
        keep_filters=args.keep_filter,
        regions=args.region,
        info_expressions=args.info_expr)


def _add_filter_args(parser):
    parser.add_argument("--min-qual", type=float, default=None,
                        help="drop records with QUAL below this value")
    parser.add_argument("--keep-filter", action="append", default=None,
                        metavar="NAME",
                        help="keep only records whose FILTER is NAME "
                             "(repeatable; 'PASS' also accepts '.')")
    parser.add_argument("--region", action="append", default=None,
                        metavar="CHR[:START[-END]]",
                        help="1-based inclusive region (repeatable)")
    parser.add_argument("--info-expr", action="append", default=None,
                        metavar="EXPR",
                        help="INFO filter, e.g. 'DP>10', 'TYPE=SNP', 'FLAG'")


def cmd_stats(args):
    with VCFReader(args.vcf) as reader:
        records = iter(reader)
        rec_filter = _build_filter(args)
        if rec_filter is not None:
            records = rec_filter.filter(records)
        collector = StatsCollector(reader.header.samples)
        for record in records:
            collector.add(record)
    print(format_report(collector))
    return 0


def cmd_filter(args):
    rec_filter = _build_filter(args)
    with VCFReader(args.vcf) as reader:
        _emit_header(reader)
        for record in rec_filter.filter(reader):
            print(record.to_line())
    return 0


def cmd_split(args):
    with VCFReader(args.vcf) as reader:
        _emit_header(reader)
        for record in split_multiallelic(reader):
            print(record.to_line())
    return 0


def _emit_header(reader):
    header = reader.header
    for key, values in header.meta.items():
        for value in values:
            print("##%s=%s" % (key, value))
    print(header.column_line)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="vcftool",
        description="Local VCF parsing, filtering, splitting, statistics. "
                    "No network access; input is local files only.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_stats = sub.add_parser("stats", help="print statistics report")
    p_stats.add_argument("vcf", help="input VCF file ('-' for stdin)")
    _add_filter_args(p_stats)
    p_stats.set_defaults(func=cmd_stats)

    p_filter = sub.add_parser("filter", help="write filtered VCF to stdout")
    p_filter.add_argument("vcf")
    _add_filter_args(p_filter)
    p_filter.set_defaults(func=cmd_filter)

    p_split = sub.add_parser("split",
                             help="split multiallelic records to stdout")
    p_split.add_argument("vcf")
    p_split.set_defaults(func=cmd_split)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    vcf = getattr(args, "vcf", None)
    if vcf == "-":
        args.vcf = sys.stdin
    try:
        return args.func(args)
    except VCFError as err:
        print("error: %s" % err, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
