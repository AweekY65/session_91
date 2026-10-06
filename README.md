# vcftool

纯本地的 VCF 变异文件解析与统计工具。所有输入（VCF、参考信息、报告）
只来自本地文件或内存流，**不访问** NCBI、Ensembl、ClinVar、云数据库或
任何外部服务。测试用 VCF 全部由测试代码现场生成。

## 支持的 VCF 范围

- VCF v4.x 文本格式（未压缩；gzip 可先解压后通过文件对象传入）
- `##` 元信息行：`fileformat`、`INFO`、`FORMAT`、`FILTER`、`contig`
- `#CHROM` 列头行与 0..n 个样本列
- 常见 SNP / INDEL 记录；多等位基因（`,` 分隔的多个 ALT）
- INFO 解析：Flag（无值）与 `key=value`；`record.info_value(key, int/float)`
  提供类型化访问
- FORMAT / sample 解析：按 `FORMAT` 键映射到每个样本的值字典
- genotype（GT）：单倍体、二倍体及更高倍性

暂不支持：结构化变异（`<DEL>` 等符号等位只做透传）、断点（BND）语义、
VCF 的索引文件（`.tbi`/`.csi`）——解析为顺序流式，不依赖索引。

## Genotype 语义

| 写法        | 含义                                   |
|-------------|----------------------------------------|
| `0/1`       | 未定相（unphased）二倍体杂合           |
| `0|1`       | 定相（phased）二倍体杂合               |
| `1/1`       | 纯合 ALT                               |
| `./.`, `.|.` | 缺失（分隔符保留定相信息）            |
| `1`         | 单倍体                                 |
| `1/2`       | 多等位记录中两个不同 ALT 的杂合        |

等位编号：`0` = REF，`1..n` = 各 ALT。判定规则：

- **called**：所有等位都非缺失
- **heterozygous**：called 且等位编号不完全相同（定相与否都算杂合）
- 统计时 called 基因型的每个非缺失等位都计入 AN/AC（单倍体计 1 个）

## 统计公式

对每个样本：

```
call rate      = called 基因型数 / 记录总数
heterozygosity = 杂合 called 基因型数 / called 基因型数
```

对每个位点（跨所有样本）：

```
AN            = 所有 called 基因型的非缺失等位总数
AC[k]         = 等位 k（0=REF, 1..n=ALT）的出现次数
AF[k]         = AC[k] / AN
```

## 过滤

`stats` 与 `filter` 子命令共享以下过滤条件（AND 组合）：

- `--min-qual Q`：QUAL ≥ Q（QUAL 为 `.` 的记录被剔除）
- `--require-pass`：FILTER 为 `PASS` 或 `.`
- `--region CHROM[:START[-END]]`：1-based 闭区间，可重复（OR）
- `--info 'KEY<op>VALUE'`：op ∈ `= == != > < >= <=`，可重复（AND）；
  两边都能解析为数字时按数值比较，否则按字符串比较

## 多等位拆分（split）

`split` 子命令把含 n 个 ALT 的记录拆成 n 条双等位记录，规则明确：

- 第 i 条输出记录只保留 ALT[i]，REF 不变
- genotype 重映射：`0 → 0`，`i → 1`，**其他等位 → 缺失 `.`**；
  倍性与定相分隔符保留
- INFO/FORMAT 按 header 声明的 `Number` 重映射：
  `Number=A` 取第 i 个值；`Number=R` 取 REF 与第 i 个值；
  `Number=G` 置为 `.`；其余原样复制

## 错误报告

所有解析错误抛出 `VCFRecordError` / `VCFHeaderError`，消息带 1-based
行号，CLI 以退出码 2 报告。覆盖：

- 列数不足（< 8 列）或与 header 样本数不符（截断行 / 多余列）
- 非法 POS（非整数或 < 1）
- FORMAT 键未在 `##FORMAT` 中声明（GT 除外）
- 样本字段值个数与 FORMAT 键数不符（截断的样本字段）
- GT 等位编号超出该记录等位范围

## 使用

```bash
# 统计（可叠加过滤条件）
python3 -m vcftool stats examples/example.vcf
python3 -m vcftool stats examples/example.vcf \
    --min-qual 50 --require-pass --region chr1:1-1000 --info 'DP>=20'

# 输出通过过滤的记录
python3 -m vcftool filter examples/example.vcf --require-pass -o out.vcf

# 拆分多等位记录
python3 -m vcftool split examples/example.vcf -o split.vcf
```

`stats` 输出三段制表符分隔结果：记录数、每样本
`call_rate/called/total/het/heterozygosity`、每位点 `AN/AC/AF`。

## 测试

```bash
python3 -m pytest tests/ -v
```

测试全部在终端运行、输出统计结果，不依赖 GUI。覆盖：SNP、INDEL、
多等位、phased/unphased、缺失 genotype、QUAL/FILTER/region/INFO 过滤、
频率统计、各类非法记录（列数、POS、未知 FORMAT、截断行）以及
5 万条记录的流式处理与无界输入的惰性读取。

## 代码结构

```
vcftool/
  parser.py    # 流式解析：header / 记录 / INFO / FORMAT / GT，带行号校验
  genotype.py  # Genotype 模型（phased、missing、多倍性）
  stats.py     # call rate、heterozygosity、AC/AN/AF
  filters.py   # QUAL / FILTER / region / INFO 表达式过滤
  split.py     # 多等位拆分与 genotype 重映射、记录序列化
  cli.py       # stats / filter / split 子命令
tests/         # pytest 测试（VCF 全部由代码生成）
examples/      # 示例 VCF
```
