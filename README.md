# TangentXY

# Overview

`TangnetXY` is a sex-informed copy number inference algorithm which extends Tangent (PMID: 36040167).

# Installation

To install `TangentXY`, you'll need either the `remotes` or `devtools` package:

```r
remotes::install_github("beroukhim-lab/RtangentXY")
```

# Usage

For more detailed function documentation, please see <https://beroukhim-lab.github.io/RtangentXY/>

The inputs to the tangent algorithm are as follows:

1. Normal sample signal matrix with values in log2(Relative Copy Number) format
2. Tumor sample signal matrix with values in log2(Relative Copy Number) format
3. Probe information file
4. Sample information file 

The package comes with example datasets to demonstrate the expected format of each input. 

## Tangent

```r
# Filter out common germline CNVs
probes_filt <- filter_blacklist(example_pif)
nsig_df <- example_nsig_df %>% filter(locus %in% probes_filt$locus)
tsig_df <- example_tsig_df %>% filter(locus %in% probes_filt$locus)

# Diagnostic plot to select number of latent factors
plot_latent_factor_importance(
  sif_df  = example_sif,
  nsig_df = nsig_df,
  tsig_df = tsig_df
)

# Run tangent
n_latent <- 5
tangent_res <- run_tangent(
  sif_df   = example_sif,
  nsig_df  = nsig_df,
  tsig_df  = tsig_df,
  n_latent = n_latent
)

# Plot signal and noise profiles before and after normalization
plot_signal_noise(
  tnorm    = tangent_res,
  sif_df   = example_sif,
  pif_df   = probes_filt,
  tsig_df  = tsig_df,
  n_latent = n_latent
)

# Plot signal profile before and after tangent normalization for a given sample
plot_signal_along_loci(
  tnorm     = tangent_res,
  sample_id = "tumor.female1",
  pif_df    = probes_filt,
  tsig_df   = tsig_df,
  n_latent  = n_latent
)
```

## Segmentation

The tangent pipeline outputs a normalized signal matrix which can then be used as input
for Circular Binary Segmentation (CBS):

```r
# Run CBS
cbs <- run_cbs(tangent_res, n_cores = 5)

# Save results to disk
save_segmentation(cbs, output_dir = ".")
```

## Pseudotangent

When relatively few normal samples are available, the pseudotangent algorithm may
be more appropriate.

```r
ptangent_res <- run_pseudotangent(
  sif_df         = example_sif,
  nsig_df        = nsig_df,
  tsig_df        = tsig_df,
  n_latent_init  = 5,
  num_partitions = 3,
  n_latent_part  = 3
)
```

# Support

Please [open an issue](https://github.com/beroukhim-lab/RtangentXY/issues/new) for support.

# Citation

These are the research papers associated with `RtangentXY`:

* Gao GF, Oh C, Saksena G, Deng D, Westlake LC, Hill BA, Reich M, Schumacher SE, Berger AC, Carter SL, Cherniack AD, Meyerson M, Tabak B, Beroukhim R, Getz G. Tangent normalization for somatic copy-number inference in cancer genome analysis. Bioinformatics. 2022 Oct 14;38(20):4677-4686. doi: 10.1093/bioinformatics/btac586. PMID: 36040167; PMCID: PMC9563697.
