"""Command line interface: stats / filter / split. Local files only."""

import argparse
import sys

from .errors import VCFError
from .filters import RecordFilter
from .parser import VCFReader
from .split import format_record, split_multiallelic
from .stats import StatsCollector


def _add_filter_args(parser):
    parser.add_argument("--min-qual", type=float, default=None,
                        help="keep records with QUAL >= this value")
    parser.add_argument("--require-pass", action="store_true",
                        help="keep only records whose FILTER is PASS or '.'")
    parser.add_argument("--region", action="append", default=[],
                        metavar="CHROM[:START[-END]]",
                        help="1-based inclusive region; repeatable")
    parser.add_argument("--info", action="append", default=[],
                        metavar="EXPR",
                        help="INFO filter expression, e.g. 'DP>=10'; "
                             "repeatable")


def _build_filter(args):
    return RecordFilter(
        min_qual=args.min_qual,
        require_pass=args.require_pass,
        regions=args.region,
        info_exprs=args.info,
    )


def _open_out(path):
    return open(path, "w") if path else sys.stdout


def _write_header(out, reader):
    for line in reader.header.meta_lines:
        out.write(line + "\n")
    out.write(reader.header.column_line + "\n")


def cmd_stats(args):
    with VCFReader(args.vcf) as reader:
        record_filter = _build_filter(args)
        collector = StatsCollector(reader.header.samples)
        for record in reader:
            if record_filter.accepts(record):
                collector.add_record(record)
    for line in collector.summary_lines():
        print(line)
    return 0


def cmd_filter(args):
    out = _open_out(args.output)
    try:
        with VCFReader(args.vcf) as reader:
            record_filter = _build_filter(args)
            _write_header(out, reader)
            kept = 0
            for record in reader:
                if record_filter.accepts(record):
                    out.write(format_record(record, reader.header.samples) + "\n")
                    kept += 1
        print("kept\t%d" % kept, file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


def cmd_split(args):
    out = _open_out(args.output)
    try:
        with VCFReader(args.vcf) as reader:
            _write_header(out, reader)
            written = 0
            for record in reader:
                if len(record.alts) <= 1:
                    out.write(format_record(record, reader.header.samples) + "\n")
                    written += 1
                else:
                    for sub in split_multiallelic(record, reader.header):
                        out.write(format_record(sub, reader.header.samples) + "\n")
                        written += 1
        print("written\t%d" % written, file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog="vcftool",
        description="Local-only VCF parser, filter, splitter and statistics.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_stats = sub.add_parser("stats", help="compute per-sample and per-site stats")
    p_stats.add_argument("vcf", help="input VCF file (local path)")
    _add_filter_args(p_stats)
    p_stats.set_defaults(func=cmd_stats)

    p_filter = sub.add_parser("filter", help="write records passing filters")
    p_filter.add_argument("vcf", help="input VCF file (local path)")
    p_filter.add_argument("-o", "--output", default=None,
                          help="output VCF path (default: stdout)")
    _add_filter_args(p_filter)
    p_filter.set_defaults(func=cmd_filter)

    p_split = sub.add_parser("split", help="split multiallelic records")
    p_split.add_argument("vcf", help="input VCF file (local path)")
    p_split.add_argument("-o", "--output", default=None,
                         help="output VCF path (default: stdout)")
    p_split.set_defaults(func=cmd_split)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except VCFError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
