"""End-to-end CLI tests: stats / filter / split print results to terminal."""

import subprocess
import sys

from conftest import HEADER

VCF_BODY = """\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\trs1\tA\tG\t50\tPASS\tDP=30\tGT\t0/1\t1/1
chr1\t200\t.\tAT\tA\t8\tq10\tDP=4\tGT\t0/1\t./.
chr1\t300\t.\tC\tT,G\t90\tPASS\tDP=40;AF=0.1,0.2\tGT\t1/2\t0/2
"""


def _write_vcf(tmp_path):
    path = tmp_path / "input.vcf"
    path.write_text(HEADER + VCF_BODY)
    return str(path)


def _run_cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "vcftool"] + list(args),
        capture_output=True, text=True, cwd=_repo_root(),
    )


def _repo_root():
    import os
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_cli_stats(tmp_path):
    result = _run_cli("stats", _write_vcf(tmp_path))
    assert result.returncode == 0
    out = result.stdout
    assert "records\t3" in out
    assert "S1\t1.0000" in out           # S1 call rate 3/3
    assert "S2\t0.6667" in out           # S2 call rate 2/3
    assert "chr1:100" in out and "chr1:300" in out


def test_cli_stats_with_filters(tmp_path):
    result = _run_cli("stats", _write_vcf(tmp_path),
                      "--min-qual", "50", "--require-pass",
                      "--region", "chr1:50-250",
                      "--info", "DP>=10")
    assert result.returncode == 0
    assert "records\t1" in result.stdout


def test_cli_filter_writes_passing_records(tmp_path):
    out_path = tmp_path / "out.vcf"
    result = _run_cli("filter", _write_vcf(tmp_path),
                      "--require-pass", "-o", str(out_path))
    assert result.returncode == 0
    lines = [l for l in out_path.read_text().splitlines()
             if not l.startswith("#")]
    assert len(lines) == 2  # q10 record removed


def test_cli_split_multiallelic(tmp_path):
    out_path = tmp_path / "split.vcf"
    result = _run_cli("split", _write_vcf(tmp_path), "-o", str(out_path))
    assert result.returncode == 0
    lines = [l for l in out_path.read_text().splitlines()
             if not l.startswith("#")]
    assert len(lines) == 4  # 3 records, one triallelic -> 2
    split_lines = [l for l in lines if l.startswith("chr1\t300\t")]
    assert len(split_lines) == 2
    assert "\tT\t" in split_lines[0] and "\tG\t" in split_lines[1]
    # genotype remapped: 1/2 -> 1/. for ALT1=T, ./1 for ALT2=G
    assert split_lines[0].endswith("\t1/.\t0/.")
    assert split_lines[1].endswith("\t./1\t0/1")


def test_cli_error_exit_code_and_line_number(tmp_path):
    bad = tmp_path / "bad.vcf"
    bad.write_text(HEADER + VCF_BODY +
                   "chr1\t400\t.\tA\tG\t50\tPASS\tDP=1\tGT\t0/1\n")
    result = _run_cli("stats", str(bad))
    assert result.returncode == 2
    assert "line" in result.stderr and "truncated" in result.stderr
