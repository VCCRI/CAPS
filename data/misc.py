"""
Functions for filtering and annotation of gnomAD variants.
See https://github.com/macarthur-lab/gnomad_lof.
"""

import hail as hl
from typing import Union


def trimer_from_heptamer(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
    return (
        t.annotate_rows(context=trimer_expr)
        if isinstance(t, hl.MatrixTable)
        else t.annotate(context=trimer_expr)
    )


def annotate_variant_types(
    t: Union[hl.MatrixTable, hl.Table], heptamers: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (
        ((t.ref == "A") & (t.alt == "G"))
        | ((t.ref == "G") & (t.alt == "A"))
        | ((t.ref == "T") & (t.alt == "C"))
        | ((t.ref == "C") & (t.alt == "T"))
    )
    cpg_expr = (
        (t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1 : mid_index] == "C")
    ) | (
        (t.ref == "C")
        & (t.alt == "T")
        & (t.context[mid_index + 1 : mid_index + 2] == "G")
    )
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (
        hl.case()
        .when(t.cpg, "CpG")
        .when(t.transition, "non-CpG transition")
        .default("transversion")
    )
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )
    else:
        return t.annotate(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return (
        hl.switch(base)
        .when("A", "T")
        .when("T", "A")
        .when("G", "C")
        .when("C", "G")
        .default(base)
    )


def reverse_complement_bases(
    bases: hl.expr.StringExpression,
) -> hl.expr.StringExpression:
    return hl.delimit(
        hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), ""
    )


def collapse_strand(
    ht: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    collapse_expr = {
        "ref": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.ref),
            ht.ref,
        ),
        "alt": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.alt),
            ht.alt,
        ),
        "context": hl.cond(
            ((ht.ref == "G") | (ht.ref == "T")),
            reverse_complement_bases(ht.context),
            ht.context,
        ),
        "was_flipped": (ht.ref == "G") | (ht.ref == "T"),
    }
    return (
        ht.annotate(**collapse_expr)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**collapse_expr)
    )


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f"[ATCG]{{{str_len}}}")
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f"[ATCG]{{{str_len}}}")
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        "methylation_level": hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }
    if annotate_coverage:
        annotation["exome_coverage"] = ht.coverage.exomes.median
    return (
        ht.annotate(**annotation)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**annotation)
    )


def get_an_adj_criteria(
    hail_table,
    sex_split,
    an_cutoff: float = 0.8,
):
    """Get lower bound allele number (AN) thresholds.

    hail_table -- variants in ht format
    sex_split -- how many males, how many females
    an_cutoff -- percent of genotype calls
    """
    return (
        hl.case()
        .when(
            hail_table.locus.in_autosome_or_par(),
            hail_table.freq[0].AN >= an_cutoff * 2 * sum(sex_split.values()),
        )
        .when(
            hail_table.locus.in_x_nonpar(),
            hail_table.freq[0].AN
            >= an_cutoff * (sex_split["male"] + sex_split["female"] * 2),
        )
        .when(
            hail_table.locus.in_y_nonpar(),
            hail_table.freq[0].AN >= an_cutoff * sex_split["male"],
        )
        .or_missing()
    )
