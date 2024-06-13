# TangentXY R Package

An easy-to-use package containing functions that run TangentXY and PseudotangentXY, sex-informed extensions of the existing Tangent and Pseudotangent inference pipelines.

## Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)
-   [Support](#support)
-   [Versions](#versions)

## Installation

Check the DESCRIPTION file for the list of all imports. Installing the latest version of `tidyverse` should take care of many of these package imports. Check [this link](https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html) for instructions on installing an R package from GitHub.

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

`filter_blacklist()` `run_tangent()` `run_cbs()` `convert_cbs()`

If the set of normals is particularly non-representative, we suggest running PseudoTangent.

`filter_blacklist()` `run_pseudptangent()`

Because CBS is an essential intermediate step within PseudoTangent, all segmentation is done within the `run_pseudotangent()` function.

## Support

Please [open an issue](https://github.com/beroukhim-lab/RtangentXY/issues/new) for support.

## Versions

0.0.1

-   Pre-release version
