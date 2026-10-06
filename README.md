# vcftool — 本地 VCF 解析与统计工具

纯本地、零网络依赖的 VCF（Variant Call Format）解析、过滤、多等位拆分与统计工具。
所有输入（VCF、参考信息、索引、报告）均来自本地文件或内存，**不会**访问 NCBI、
Ensembl、ClinVar、云数据库或任何外部服务。测试所用 VCF 全部由工程自带
（`tests/data/`）或由测试代码在内存中生成。

## 环境要求

- Python 3.8+（仅使用标准库，无第三方运行时依赖）
- 运行测试需要 `pytest`

## 支持的 VCF 范围

- VCF v4.2 语法的常见子集：`##` 元信息行（含 `INFO` / `FORMAT` / `contig`
  结构化声明）、`#CHROM` 列头行、数据行。
- 变异类型：SNP、INDEL（插入/缺失）、多等位（multiallelic）记录；支持符号
  ALT 之外的常规 ACGTN 等位基因。
- INFO 字段按 header 中声明的 `Number`（`0/1/A/R/G/.`）与 `Type`
  （`Integer/Float/Flag/Character/String`）解析；`Number=A/R/G` 统一解析为列表。
- FORMAT 字段同样按 header 声明解析；`GT` 被隐式识别，其余未在 header 中
  声明的 FORMAT 键会报“未知 FORMAT”错误。
- 样本列为可选（支持无样本的 sites-only VCF）。
- 流式（streaming）解析：逐行读取，内存占用与记录数无关，可处理大文件。

不支持的特性：结构化/断点变异（BND）、符号 ALT 的语义解释、VCF 压缩
（bgzip/.gz 需先解压）、参考基因组一致性校验（不读取任何外部参考文件）。

## Genotype 语义

- `0/1`、`1|0` 等：等位基因索引列表，`0` 为 REF，`1..n` 为各 ALT。
- `|` 表示 phased，`/` 表示 unphased；同一 GT 中混用分隔符按 unphased 处理。
- `.`、`./.`、`.|.` 为缺失 genotype；`0/.` 这类部分缺失仍视为已检出（called），
  缺失的等位基因不计入 AN/AC。
- 倍性（ploidy）由等位基因个数决定，支持单倍体与二倍体及更高倍性。
- 杂合（heterozygous）判定：已检出的等位基因中存在两个及以上不同取值。

## 统计公式

按样本（对整个记录流累计）：

- `call_rate = called_genotypes / total_records`
- `heterozygosity = heterozygous_calls / called_genotypes`（仅统计二倍体杂合）

按记录（等位基因索引 0 = REF，1..n = 各 ALT）：

- `AC`（allele count）：所有样本 GT 中等于该索引的非缺失等位基因个数
- `AN`（allele number）：所有样本 GT 中非缺失等位基因总数
- `AF`（allele frequency）：`AC / AN`；`AN = 0` 时为 `NA`

## 多等位拆分规则（split）

每个 ALT 等位基因输出一条独立的二倍等位记录，规则如下：

- `GT`：等位 `0` → `0`；等于当前 ALT（索引 i+1）→ `1`；其余 ALT → 缺失 `.`；
  缺失保持缺失；倍性与 phased 状态保留。
- INFO/FORMAT `Number=A`：取第 i 个元素；`Number=R`：取第 0 与第 i+1 个元素；
  `Number=G`：丢弃（genotype 似然无法安全重映射）；`Number=1/.`：原样保留。
- `CHROM/POS/ID/REF/QUAL/FILTER` 原样复制。

## 过滤

支持任意组合的过滤条件：

- `--min-qual Q`：QUAL 低于 Q（或 QUAL 缺失为 `.`）的记录被丢弃。
- `--keep-filter NAME`：仅保留 FILTER 列等于 NAME 的记录（可重复；`PASS`
  同时接受 `.`）。
- `--region CHR[:START[-END]]`：1-based 闭区间，按变异跨度（POS..POS+len(REF)-1）
  与区间重叠判定（可重复，并集）。
- `--info-expr EXPR`：INFO 表达式，如 `DP>10`、`AF<=0.05`、`TYPE=SNP`、
  `TYPE!=INDEL`、裸键 `DB`（存在性）；列表值按“任一元素满足”判定。

## 命令行用法

```bash
# 统计报告（可叠加过滤条件）
python3 -m vcftool stats input.vcf --min-qual 30 --region chr1:1000-2000 \
    --keep-filter PASS --info-expr "DP>10"

# 输出过滤后的 VCF 到 stdout
python3 -m vcftool filter input.vcf --region chr2

# 拆分多等位记录到 stdout
python3 -m vcftool split input.vcf

# 从 stdin 读取
cat input.vcf | python3 -m vcftool stats -
```

解析错误（列数错误、非法 POS、未知 FORMAT、截断行等）会以非零退出码终止，
并在 stderr 输出带行号的错误信息，例如：
`error: line 3: invalid POS 'not_a_pos': not an integer`。

## Python API

```python
from vcftool import VCFReader, RecordFilter, collect_stats, split_multiallelic

with VCFReader("input.vcf") as reader:          # 也可传入字符串行迭代器
    records = RecordFilter(min_qual=30).filter(reader)
    for rec in split_multiallelic(records):
        print(rec.chrom, rec.pos, rec.alts, rec.samples[0]["GT"])
```

## 运行测试

```bash
python3 -m pytest tests/ -v
```

测试全部在终端执行并输出统计结果，不依赖 GUI，覆盖：SNP、INDEL、多等位记录、
phased/unphased 与缺失 genotype、各类过滤条件、等位基因频率统计、非法记录
（列数错误 / 非法 POS / 未知 FORMAT / 截断行，均校验行号）、多等位拆分与
genotype 重映射，以及 2 万条记录的流式处理（含内存占用上界断言）。

## 目录结构

```
vcftool/
  parser.py     # 流式 VCF 解析器（header / INFO / FORMAT / GT）
  genotype.py   # GT 语义：phased、缺失、杂合判定
  stats.py      # call rate、杂合度、AC/AN/AF 统计
  filters.py    # QUAL / FILTER / 区间 / INFO 表达式过滤
  split.py      # 多等位拆分与 genotype 重映射
  cli.py        # 命令行入口（stats / filter / split）
tests/
  data/         # 工程自带的测试 VCF
  test_*.py     # 各场景自动化测试
```
