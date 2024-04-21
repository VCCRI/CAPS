# Context-Adjusted Proportion of Singletons (CAPS)

## Code for the paper

### The "Data" pipeline

`data/` contains [Hail](https://hail.is/) and Snakemake code that requires execution in Google Cloud and saves files to a Google Storage bucket. Copies of the generated files are available in `files/`.

#### How to run

1. Create a new cluster: `hailctl dataproc start <cluster_name> --packages snakemake --requester-pays-allow-buckets gnomad-public-requester-pays --project <project_name> --bucket <bucket_name> --region <region> --num-workers <N> --image-version=2.0.27-debian10`
2. Connect to the cluster: `gcloud beta compute ssh <user_name>@<cluster_name>-m --zone "<zone>" --project "<project_name>"`
3. `git clone` this repository and navigate to `data/`
4. Run the pipeline: `snakemake --cores all --configfile config.yaml --config gcp_rootdir="<bucket_name>/some_directory/" gcp_username="<user_name>"`

### The "Analysis" pipeline

`analysis/` contains scripts that calculate and visualise CAPS scores using files created in `data/`.

#### How to run

1. Navigate to `CAPS/analysis/`
2. `snakemake --cores all --config gcp="False"` (faster: uses copies from `files/`) or `snakemake --cores all --config gcp="True" gcp_rootdir="<bucket_name>/some_directory/"` (slower: uses GS files)

## Using CAPS

### Custom sets of variants

To get CAPS estimates for your set of variants, use the `template` file: `snakemake -s template -c1 -C [KEY=VALUE ...]`. The required values are
- `obs` (grouped variants annotated with at least `context`, `ref`, `alt`, `methylation_level`, `singleton_count` and `variant_count` fields)
- `exp` (expected proportions, one per `context`-`ref`-`alt`-`methylation_level` group)
- `var` (variable of interest, must be a valid field in `obs`)
- `calculate_caps_script` (`calculate_caps.R`)
- `viz_scores_script` (`viz_scores.R`)
- `scores` (filename for the output scores)
- `plot` (filename for the output plot)

For example, `snakemake -s template -c1 -C obs=analysis/canonical_splice_site_vars.tsv exp=model/phat.tsv var=worst_csq calculate_caps_script=analysis/calculate_caps.R viz_scores_script=analysis/viz_scores.R scores=scores.tsv plot=plot.pdf`.
