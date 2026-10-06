"""Split multiallelic records into one biallelic record per ALT allele.

Genotype remapping rule (per output record for ALT index i, 1-based):
  * allele 0 (REF)  -> 0
  * allele i        -> 1
  * any other allele -> missing ('.')
Ploidy and phasing are preserved.

INFO/FORMAT value remapping uses the Number declared in the header:
  * Number=A: keep only the value for ALT i
  * Number=R: keep the REF value and the value for ALT i
  * Number=G: dropped ('.') -- genotype-order values are not meaningful
  * anything else (scalar/undeclared): copied unchanged
"""

from .parser import VCFRecord


def _remap_list(value, indices):
    parts = value.split(",")
    out = []
    for idx in indices:
        if idx >= len(parts):
            return None
        out.append(parts[idx])
    return ",".join(out)


def _remap_numbered(value, number, alt_idx):
    if value is True:
        return value
    if number == "A":
        return _remap_list(value, [alt_idx])
    if number == "R":
        return _remap_list(value, [0, alt_idx + 1])
    if number == "G":
        return None
    return value


def _remap_gt(gt_text, alt_idx):
    from .genotype import Genotype
    gt = Genotype.parse(gt_text)
    new_alleles = []
    for allele in gt.alleles:
        if allele is None:
            new_alleles.append(None)
        elif allele == 0:
            new_alleles.append(0)
        elif allele == alt_idx + 1:
            new_alleles.append(1)
        else:
            new_alleles.append(None)
    return str(Genotype(new_alleles, gt.phased))


def split_multiallelic(record, header):
    """Yield one biallelic VCFRecord per ALT allele of `record`."""
    for alt_idx, alt in enumerate(record.alts):
        info = {}
        for key, value in record.info.items():
            decl = header.info.get(key)
            number = decl.get("Number") if decl else None
            remapped = _remap_numbered(value, number, alt_idx)
            if remapped is not None:
                info[key] = remapped

        sample_data = []
        for data in record.sample_data:
            new_data = {}
            for key, value in data.items():
                if key == "GT":
                    new_data[key] = _remap_gt(value, alt_idx)
                    continue
                decl = header.formats.get(key)
                number = decl.get("Number") if decl else None
                remapped = _remap_numbered(value, number, alt_idx)
                new_data[key] = remapped if remapped is not None else "."
            sample_data.append(new_data)

        yield VCFRecord(
            record.chrom, record.pos, record.vid, record.ref, [alt],
            record.qual, list(record.filters), info,
            list(record.format_keys), sample_data, record.line_no,
        )


def format_record(record, sample_names):
    """Render a record back to a VCF text line (for filter/split output)."""
    cols = [
        record.chrom, str(record.pos), record.vid, record.ref,
        ",".join(record.alts) if record.alts else ".",
        "." if record.qual is None else ("%g" % record.qual),
        ";".join(record.filters) if record.filters else ".",
    ]
    if record.info:
        items = []
        for key, value in record.info.items():
            items.append(key if value is True else "%s=%s" % (key, value))
        cols.append(";".join(items))
    else:
        cols.append(".")
    if sample_names:
        cols.append(":".join(record.format_keys))
        for data in record.sample_data:
            cols.append(":".join(data.get(k, ".") for k in record.format_keys))
    return "\t".join(cols)
