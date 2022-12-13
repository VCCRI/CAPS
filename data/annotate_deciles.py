import hail as hl
def annotate_deciles(ht, variable):
    """Group into decile bins.

    ht -- variants (Hail table)
    variable -- variable of interest
    """

    quantiles = ht.aggregate(
        hl.agg.approx_quantiles(
            ht[variable], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        )
    )
    ht = ht.annotate(
        variable_bin=hl.case()
        .when(ht[variable] >= quantiles[9], 9)
        .when(ht[variable] >= quantiles[8], 8)
        .when(ht[variable] >= quantiles[7], 7)
        .when(ht[variable] >= quantiles[6], 6)
        .when(ht[variable] >= quantiles[5], 5)
        .when(ht[variable] >= quantiles[4], 4)
        .when(ht[variable] >= quantiles[3], 3)
        .when(ht[variable] >= quantiles[2], 2)
        .when(ht[variable] >= quantiles[1], 1)
        .when(ht[variable] >= quantiles[0], 0)
        .or_missing()
    )
    return ht
