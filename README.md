# TangentXY R Package

An easy-to-use package containing functions that run TangentXY and PseudotangentXY, sex-informed extensions of the existing Tangent and Pseudotangent inference pipelines.

## Table of Contents

-   [Cite](#cite)
-   [Installation](#installation)
-   [Usage](#usage)
-   [Support](#support)
-   [Versions](#versions)


## Cite

These are the research papers associated with `RtangentXY`:

* Gao GF, Oh C, Saksena G, Deng D, Westlake LC, Hill BA, Reich M, Schumacher SE, Berger AC, Carter SL, Cherniack AD, Meyerson M, Tabak B, Beroukhim R, Getz G. Tangent normalization for somatic copy-number inference in cancer genome analysis. Bioinformatics. 2022 Oct 14;38(20):4677-4686. doi: 10.1093/bioinformatics/btac586. PMID: 36040167; PMCID: PMC9563697.

## Installation

Check the `RtangentXY/DESCRIPTION` file for the list of all imports. This is an R package, so it should be able to run on any device with R installed. We recommend using the latest version of R, as some of the dependencies may require more recent versions of R. Check [this link](https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html) for instructions on installing an R package from GitHub.

We will present one way to install this package. First download this repository into your workspace:

```sh
git clone git@github.com:beroukhim-lab/RtangentXY.git
cd RtangentXY
```

You can then install `RtangentXY` by running the following commands in an R environment, within the `/RtangentXY` directory:

```r
install.packages("devtools")
devtools::install()
```

You can then test if the package and dependencies were properly installed by running the following command in an R environment, within the `/RtangentXY` directory:

```r
devtools::test()
```

## Expected Inputs

SIF File

PIF File

Tumor Signal File

Normal Signal File

## Usage

Once this package is installed, load the package by calling:

`library(RtangentXY)`

We suggest the following workflows:

If you have a subset of normal samples, we suggest running Tangent.

`filter_blacklist()` `run_tangent()` `run_cbs()` 

If the set of normals is particularly non-representative, we suggest running PseudoTangent.

`filter_blacklist()` `run_pseudptangent()`

Because CBS is an essential intermediate step within PseudoTangent, all segmentation is done within the `run_pseudotangent()` function.

Both of these workflows will output a DNAcopy object, which can be plotted.

## Support

Please [open an issue](https://github.com/beroukhim-lab/RtangentXY/issues/new) for support.

## Versions

0.0.1

-   Pre-release version
