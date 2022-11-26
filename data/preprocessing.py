import hail as hl

from misc import *

def preprocessing(
    exomes_ht,
    context_ht,
    sex_split,
    msc_type,
    biotype="protein_coding",
):
    """Preprocessing steps for selected variants.

    exomes_ht -- WES variants (Hail table)
    context_ht -- context (Hail table)
    sex_split -- how many males, how many females
    msc_type -- variant type based on most severe consequence
    biotype -- which transcripts should be selected
    """

    exomes = hl.read_table(exomes_ht)

    # Allele number (AN) adjustment.
    exomes = exomes.filter(get_an_adj_criteria(exomes, sex_split))

    # Filter the table so that only those variants that have AF>0 and
    # filter PASS are retained.  The first condition is necessary
    # because in gnomAD variants that were excluded from the analysis
    # through QC have AF=0. The condition on "most_severe_consequence"
    # removes all pLoF variants, including variants in canonical
    # splice sites.
    variants = exomes.filter(
        (exomes.freq[0].AF > 0)
        & (exomes.filters.length() == 0)
        & (exomes.vep.variant_class == "SNV")
        & (exomes.vep.most_severe_consequence == msc_type)
    )

    # For each variant there is a list of
    # "transcript_consequences". The code below extracts information
    # about the first transcript matching "biotype" where the variant
    # has the type "msc_type" (there may be several such
    # transcripts). Those mutations with "msc_type" as their most
    # severe consequence that don't have such a transcript will become
    # NAs.
    variants = variants.annotate(
        transcript_consequences=variants.vep.transcript_consequences.find(
            lambda x: (x.consequence_terms == [msc_type]) & (x.biotype == biotype)
        )
    )

    # Methylation and other context data.
    context = hl.read_table(context_ht)
    context = context[variants.key]
    # The 2020 version of MAPS uses methylation.
    # Function "prepare_ht" annotates the input table with methylation level,
    # coverage (optional), CpG/Non-CpG info, context for mutability
    # (ref allele in the middle plus two bases to the left and to the right)
    # and other, less important information.
    # For example, a variant that has changed "|...|..t|Cta|...|"
    # to "|...|..t|Tta|...|" will have "ref_codon"="Cta",
    # "alt_codon"="Tta", "ref"="C", "alt"="T" and "context"="TCT".
    variants = prepare_ht(
        variants.annotate(context=context.context, methylation=context.methylation),
        trimer=True,
        annotate_coverage=False,
    )

    return variants
