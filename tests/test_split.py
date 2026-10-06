"""Multiallelic splitting and genotype remapping."""

from vcftool import VCFReader, split_multiallelic


def _split_one(vcf_factory, record_line):
    reader = VCFReader(vcf_factory([record_line]))
    record = next(reader)
    return list(split_multiallelic(record, reader.header))


def test_split_biallelic_passthrough(vcf_factory):
    rec = "chr1\t1\t.\tA\tG\t10\tPASS\tDP=5\tGT\t0/1\t0/0"
    reader = VCFReader(vcf_factory([rec]))
    record = next(reader)
    subs = list(split_multiallelic(record, reader.header))
    assert len(subs) == 1
    assert subs[0].alts == ["G"]


def test_split_triallelic_genotype_remap(vcf_factory):
    rec = "chr1\t100\t.\tA\tG,T\t60\tPASS\tDP=9\tGT\t0/1\t1/2"
    subs = _split_one(vcf_factory, rec)
    assert len(subs) == 2
    assert [s.alts for s in subs] == [["G"], ["T"]]
    # ALT1=G record: 0/1 stays 0/1; allele 2 becomes missing -> 1/.
    gts1 = subs[0].genotypes()
    assert str(gts1[0]) == "0/1"
    assert str(gts1[1]) == "1/."
    # ALT2=T record: allele 1 becomes missing -> 0/.; 1/2 -> ./1
    gts2 = subs[1].genotypes()
    assert str(gts2[0]) == "0/."
    assert str(gts2[1]) == "./1"


def test_split_preserves_phasing(vcf_factory):
    rec = "chr1\t100\t.\tA\tG,T\t60\tPASS\tDP=9\tGT\t0|2\t1|1"
    subs = _split_one(vcf_factory, rec)
    # ALT1=G: 0|2 -> 0|. (allele 2 missing), 1|1 -> 1|1
    assert str(subs[0].genotypes()[0]) == "0|."
    assert str(subs[0].genotypes()[1]) == "1|1"
    # ALT2=T: 0|2 -> 0|1, 1|1 -> .|.
    assert str(subs[1].genotypes()[0]) == "0|1"
    assert str(subs[1].genotypes()[1]) == ".|."


def test_split_info_number_a_and_r(vcf_factory):
    rec = ("chr1\t100\t.\tA\tG,T\t60\tPASS\tDP=9;AF=0.1,0.2\t"
           "GT:AD\t0/1:5,3,4\t0/2:5,3,4")
    subs = _split_one(vcf_factory, rec)
    # AF is Number=A: each split keeps its own value; DP scalar kept.
    assert subs[0].info["AF"] == "0.1"
    assert subs[1].info["AF"] == "0.2"
    assert subs[0].info["DP"] == "9"
    # AD is Number=R: keep REF depth + the ALT's depth.
    assert subs[0].sample_data[0]["AD"] == "5,3"  # REF + ALT1 depths
    assert subs[1].sample_data[0]["AD"] == "5,4"  # REF + ALT2 depths
    assert subs[1].sample_data[1]["AD"] == "5,4"
