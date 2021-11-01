
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VariantsID

<!-- badges: start -->
<!-- badges: end -->

The goal of VariantsID is to …

## Installation

You can install the development version of VariantsID like so:

``` r
devtools::install_github("Linda24bc/VariantsID")
```

## Introduction

-   Identify the Hb variants by comparing the experimental data in the
    database
-   Original database contains the theoretical diagnostic ions of Hb
    variants been tested by mass spectrometry so far

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(VariantsID)
```

### Step 1. Input the database including the diagnostic ions of Hb varints and use MS1 data to narrow down the database - subset the database

1.1 Input the original database

HbDatabase &lt;- read\_csv(“Hb Variants\_OriginalDatabase.csv”)

1.2 Use the MS1 data to narrow down the database, if the mass shift is
about -0.93 Da, then the Mshift is -0.93 Da and the error tolerence is
0.06 Da. Thoese two values are changable and depend on the accuracy of
deconvolution.

ref &lt;- SubDatabase(HbDatabase, Mshift= -0.93, error\_Da\_L=-0.05,
error\_Da\_R=0.06)

### Step 2. Input the deconvolve MS2 results

The list should contain two columns, Exp\_m/z vs Exp\_Intensity)

exp &lt;- read\_csv(“expt mass\_AC.csv”)

### Step 3. Search the experimental results in the subset database with Variant Identifier

Run the function Variants.Identifier, the ppm\_error range is changable
and depends on the accuracy of the MS2 data. View the result list and
get the identification.

ID.results &lt;- Variants.Identifier(ref, exp, ppm\_error\_start=-2,
ppm\_error\_end=5)

### Step 4. Output the results in .csv

write.csv(ID.results, “ID\_HbAE\_1.csv”, row.names = FALSE)
