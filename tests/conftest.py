import os

import pytest

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def data_path():
    def _get(name):
        return os.path.join(DATA_DIR, name)
    return _get


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""


def make_vcf_lines(records):
    """Build an in-memory VCF (list of lines) from record column strings."""
    return (HEADER + "\n".join(records) + "\n").splitlines(True)


@pytest.fixture
def vcf_lines():
    return make_vcf_lines
