import subprocess
import sys


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "vcftool"] + list(args),
                          capture_output=True, text=True)


def test_cli_stats(data_path, capsys=None):
    result = run_cli("stats", data_path("basic.vcf"))
    assert result.returncode == 0, result.stderr
    assert "== Variant summary ==" in result.stdout
    assert "records:        5" in result.stdout
    assert "HG001" in result.stdout


def test_cli_stats_with_filters(data_path):
    result = run_cli("stats", data_path("basic.vcf"),
                     "--min-qual", "40", "--region", "chr1:1-500")
    assert result.returncode == 0, result.stderr
    assert "records:        2" in result.stdout


def test_cli_filter(data_path):
    result = run_cli("filter", data_path("basic.vcf"), "--keep-filter", "PASS")
    assert result.returncode == 0, result.stderr
    body = [l for l in result.stdout.splitlines() if not l.startswith("#")]
    assert len(body) == 3


def test_cli_split(data_path):
    result = run_cli("split", data_path("basic.vcf"))
    assert result.returncode == 0, result.stderr
    body = [l for l in result.stdout.splitlines() if not l.startswith("#")]
    assert len(body) == 6
    assert all(len(l.split("\t")[4].split(",")) == 1 for l in body)


def test_cli_error_exit_code(data_path, tmp_path):
    bad = tmp_path / "bad.vcf"
    bad.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\tnot_a_pos\t.\tA\tG\t10\tPASS\t.\n")
    result = run_cli("stats", str(bad))
    assert result.returncode == 2
    assert "line 3" in result.stderr
    assert "invalid POS" in result.stderr
