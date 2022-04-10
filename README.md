
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VariantsID

<!-- badges: start -->
<!-- badges: end -->

The goal of VariantsID is to identify Hb variants by deconvoluved MS2 data and predict diagnostic product ions of Hb beta variants by referring to the experimentally determined fragments of HBA beta

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

#### library(VariantsID)

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

exp &lt;- read\_csv(“expt mass\_cHbSS.csv”)

### Step 3. Search the experimental results in the subset database with Variant Identifier

Run the function Variants.Identifier, the ppm\_error range is changable
and depends on the accuracy of the MS2 data. View the result list and
get the identification.

ID.results &lt;- Variants.Identifier(ref, exp, ppm\_error\_start=-2,
ppm\_error\_end=5)

### Step 4. Output the results in .csv

write.csv(ID.results, “ID\_cHbSS.csv”, row.names = FALSE)

# PredictDiag

## Introduction

-   This procedure can predict the diagnostic ions of Hb beta variants
    that have been experimentally teseted by refering to the product
    ions of Hb A beta.
-   The output list from step 5 can be combined with the original
    databse in Variants Identifier as the updated database for
    sereaching.

### Step 1: Input the list of residue numbers of possible diagnostic ions for each AA in the Hb beta,and the list of reference product ions for HbA beta

diag\_ref &lt;- read.csv(“finddiag.csv”)

WT\_ref &lt;- read.csv(“ref mass list\_pro\_1.csv”)

### Step 2: Input the sequennces of HbA beta and Hb beta variants

Multiple sequences of variants sequences can be included in one .fasta
file, the sequences should have the N-terminal Met while the comparison
results exclude the N-ternimal Met.

Hbvariants &lt;- seqinr::read.fasta(file = “Hbvariants.fasta”, seqtype =
“AA”,as.string = FALSE)

WT &lt;- seqinr::read.fasta(file = “HbA.fasta”, seqtype = “AA”,as.string
= FALSE)

### Step 3: Predict the diagnostic ions by running the function

PD.result &lt;- PredictDiag(WT,WT\_ref,diag\_ref,Hbvarinats)

### Step 4:Output results in .csv file

write.csv(PD.result, “PredictDiag\_variants20.csv”, row.names = FALSE)
