## code to prepare `testData` dataset goes here


library(readr)

ref <- read.csv("HbVariants_OriginalDatabase.csv")
exp <- read.csv("expt mass_cHbSS.csv")


diag_ref <- read.csv("finddiag.csv")
WT_ref <- read.csv("ref mass list_pro_1.csv")
Hbvariants <- seqinr::read.fasta(file = "Hbvariants.fasta", seqtype = "AA",as.string = FALSE)
WT <- seqinr::read.fasta(file = "HbA.fasta", seqtype = "AA",as.string = FALSE)

