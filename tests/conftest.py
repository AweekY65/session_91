import io

import pytest

HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DB,Number=0,Type=Flag,Description="db membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
"""


def make_vcf(records, samples=("S1", "S2"), header=HEADER):
    """Build an in-memory VCF text stream from record strings."""
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        cols += "\tFORMAT\t" + "\t".join(samples)
    text = header + cols + "\n" + "\n".join(records)
    if records:
        text += "\n"
    return io.StringIO(text)


@pytest.fixture
def vcf_factory():
    return make_vcf
