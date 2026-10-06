from vcftool import VCFReader, split_multiallelic
from vcftool.genotype import Genotype


def load_split(path):
    with VCFReader(path) as reader:
        return list(split_multiallelic(reader))


def test_split_count_and_biallelic(data_path):
    records = load_split(data_path("basic.vcf"))
    # original 5 records, one triallelic -> 6 records
    assert len(records) == 6
    assert all(not r.is_multiallelic for r in records)


def test_split_gt_remapping(data_path):
    records = load_split(data_path("basic.vcf"))
    sub1, sub2 = records[2], records[3]  # from chr1:300 C>T,G
    assert sub1.alts == ["T"] and sub2.alts == ["G"]
    # HG002 was 1/2: alt1 (T) -> 1/. ; alt2 (G) -> ./1
    assert sub1.samples[1]["GT"] == Genotype([1, None])
    assert sub2.samples[1]["GT"] == Genotype([None, 1])
    # HG003 was 2/2: alt1 -> ./. ; alt2 -> 1/1
    assert sub1.samples[2]["GT"].is_missing
    assert sub2.samples[2]["GT"] == Genotype([1, 1])


def test_split_preserves_phase(data_path):
    with VCFReader(data_path("phased_multi.vcf")) as reader:
        subs = list(split_multiallelic(reader))
    assert str(subs[0].samples[0]["GT"]) == "0|."
    assert str(subs[1].samples[0]["GT"]) == "0|1"


def test_split_info_number_a(data_path):
    records = load_split(data_path("basic.vcf"))
    sub1, sub2 = records[2], records[3]
    assert sub1.info["AF"] == 0.3
    assert sub2.info["AF"] == 0.1
    assert sub1.info["TYPE"] == "SNP"
    assert sub1.info["DB"] is True
    assert sub1.info["DP"] == 40  # Number=1 kept


def test_split_format_number_r(data_path):
    records = load_split(data_path("basic.vcf"))
    sub1, sub2 = records[2], records[3]
    # AD was Number=R: 10,5,0 -> alt1: 10,5 ; alt2: 10,0
    assert sub1.samples[0]["AD"] == [10, 5]
    assert sub2.samples[0]["AD"] == [10, 0]
    assert sub1.samples[1]["AD"] == [8, 6]
    assert sub2.samples[1]["AD"] == [8, 6]


def test_split_drops_number_g(data_path):
    with VCFReader(data_path("gl_multi.vcf")) as reader:
        subs = list(split_multiallelic(reader))
    assert "PL" not in subs[0].samples[0]
    assert "PL" not in subs[0].format_keys
    assert subs[0].samples[0]["GT"] is not None


def test_split_roundtrip_line(data_path):
    records = load_split(data_path("basic.vcf"))
    line = records[2].to_line()
    cols = line.split("\t")
    assert cols[4] == "T"
    assert "AF=0.3" in cols[7]
