"""Split multiallelic records into one biallelic record per ALT allele.

Remapping rules (applied for ALT index i, 0-based):
  GT:      allele 0 -> 0; allele i+1 -> 1; any other ALT -> missing (.);
           missing stays missing. Ploidy and phasing are preserved.
  INFO/FORMAT Number=A: keep element i (ALT-ordered).
  INFO/FORMAT Number=R: keep elements [0, i+1] (REF + this ALT).
  Number=G: dropped (genotype likelihoods cannot be remapped safely).
  Number=1 or '.': kept unchanged.
QUAL, FILTER, ID, CHROM, POS, REF are copied unchanged.
"""

from .genotype import Genotype
from .parser import VariantRecord


def _remap_gt(gt, alt_index):
    target = alt_index + 1
    alleles = []
    for allele in gt.alleles:
        if allele is None:
            alleles.append(None)
        elif allele == 0:
            alleles.append(0)
        elif allele == target:
            alleles.append(1)
        else:
            alleles.append(None)
    return Genotype(alleles, phased=gt.phased)


def _subset_value(value, fdef, alt_index):
    """Subset a Number=A/R/G value for one ALT; returns _DROP to remove."""
    if fdef is None or not isinstance(value, list):
        return value
    if fdef.number == "A":
        return value[alt_index] if alt_index < len(value) else None
    if fdef.number == "R":
        picked = []
        for idx in (0, alt_index + 1):
            picked.append(value[idx] if idx < len(value) else None)
        return picked
    if fdef.number == "G":
        return _DROP
    return value


class _Drop(object):
    pass


_DROP = _Drop()


def split_record(record):
    """Yield one biallelic VariantRecord per ALT allele."""
    if not record.is_multiallelic:
        yield record
        return
    header = record.header
    for alt_index, alt in enumerate(record.alts):
        info = {}
        for key, value in record.info.items():
            new_value = _subset_value(value, header.info_defs.get(key),
                                      alt_index)
            if new_value is not _DROP:
                info[key] = new_value
        samples = []
        for sample in record.samples:
            new_sample = {}
            for key, value in sample.items():
                if key == "GT":
                    new_sample[key] = _remap_gt(value, alt_index)
                else:
                    new_value = _subset_value(
                        value, header.format_def(key), alt_index)
                    if new_value is not _DROP:
                        new_sample[key] = new_value
            samples.append(new_sample)
        format_keys = [k for k in record.format_keys
                       if k == "GT" or any(k in s for s in samples)]
        yield VariantRecord(
            record.chrom, record.pos, record.id, record.ref, [alt],
            record.qual, record.filter, info, format_keys, samples,
            record.line_number, header)


def split_multiallelic(records):
    for record in records:
        for sub in split_record(record):
            yield sub
