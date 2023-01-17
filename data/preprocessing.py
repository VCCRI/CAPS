import hail as hl

# from gnomad_lof.constraint_utils.generic import get_an_adj_criteria
# from gnomad_lof.constraint_utils.constraint_basics import prepare_ht
# from gnomad_hail.utils.generic import filter_vep_to_canonical_transcripts
# from gnomad_lof.constraint.summary_statistics import (
#     get_worst_consequence_with_non_coding,
# )

from misc import (
    get_an_adj_criteria,
    prepare_ht,
    filter_vep_to_canonical_transcripts,
    get_worst_consequence_with_non_coding,
)


def preprocessing(
    data_ht,
    context_ht,
    mutation_rates_ht,
    coverage_ht,
    sex_split,
):
    """Preprocessing steps for selected variants.

    data_ht -- WES or WGS variants (Hail table)
    context_ht -- context (Hail table)
    mutation_rates_ht -- mutability (Hail table)
    coverage_ht -- coverage (Hail table)
    sex_split -- how many males, how many females
    """

    context_ht = hl.read_table(context_ht)
    mutation_rates_ht = hl.read_table(mutation_rates_ht)
    coverage_ht = hl.read_table(coverage_ht)

    ht = hl.read_table(data_ht)

    ht = ht.annotate(coverage=coverage_ht[ht.locus].median)

    ht = ht.filter((hl.len(ht.filters) == 0) & get_an_adj_criteria(ht, sex_split))

    ht = filter_vep_to_canonical_transcripts(ht)
    ht = get_worst_consequence_with_non_coding(ht)

    context = context_ht[ht.key]
    ht = prepare_ht(
        ht.annotate(context=context.context, methylation=context.methylation),
        True,
        False,
    )

    ht = ht.annotate(
        mu=mutation_rates_ht[
            hl.struct(
                context=ht.context,
                ref=ht.ref,
                alt=ht.alt,
                methylation_level=ht.methylation_level,
            )
        ].mu_snp
    )

    return ht
