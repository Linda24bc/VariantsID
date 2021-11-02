#' PredictDiag
#'
#' @param WT sequence of wild-type(HbA beta) protein .fasta file
#' @param WT_ref reference list of fragments for wild-type protein (HbA beta)
#' @param diag_ref possible diagnostic ions for each AA of Hba beta
#' @param Hbvarinats sequences of Hb variants .fasta file
#' @importFrom   dplyr filter left_join select between
#' @importFrom   data.table as.data.table
#' @importFrom   tidyr separate
#' @importFrom  seqinr getSequence
#' @importFrom magrittr %>%



PredictDiag <- function(WT,WT_ref,diag_ref,Hbvarinats) {
  #monomz function
  monomz <- function (sequence, fragments = "by")
  {
    results_list <- vector("list")
    for (sequence_number in 1:length(sequence)) {
      peptide_vector <- strsplit(sequence[sequence_number], split = "")[[1]]
      peptide_length <- length(peptide_vector)
      if (peptide_length < 2)
        stop("sequence must contain two or more residues")
      C <- 12
      H <- 1.007825035
      O <- 15.99491463
      S <- 31.9720707
      P <- 30.973762
      N <- 14.0030740
      proton <- 1.0072764668
      residueMass <- function(residue) {
        if (residue == "A")
          mass = C * 3 + H * 5 + N + O
        if (residue == "R")
          mass = C * 6 + H * 12 + N * 4 + O
        if (residue == "N")
          mass = C * 4 + H * 6 + N * 2 + O * 2
        if (residue == "D")
          mass = C * 4 + H * 5 + N + O * 3
        if (residue == "E")
          mass = C * 5 + H * 7 + N + O * 3
        if (residue == "Q")
          mass = C * 5 + H * 8 + N * 2 + O * 2
        if (residue == "G")
          mass = C * 2 + H * 3 + N + O
        if (residue == "H")
          mass = C * 6 + H * 7 + N * 3 + O
        if (residue == "I")
          mass = C * 6 + H * 11 + N + O
        if (residue == "L")
          mass = C * 6 + H * 11 + N + O
        if (residue == "K")
          mass = C * 6 + H * 12 + N * 2 + O
        if (residue == "M")
          mass = C * 5 + H * 9 + N + O + S
        if (residue == "F")
          mass = C * 9 + H * 9 + N + O
        if (residue == "P")
          mass = C * 5 + H * 7 + N + O
        if (residue == "S")
          mass = C * 3 + H * 5 + N + O * 2
        if (residue == "T")
          mass = C * 4 + H * 7 + N + O * 2
        if (residue == "W")
          mass = C * 11 + H * 10 + N * 2 + O
        if (residue == "Y")
          mass = C * 9 + H * 9 + N + O * 2
        if (residue == "V")
          mass = C * 5 + H * 9 + N + O
        if (residue == "C")
          mass = C * 3 + H * 5 + N + O + S

        return(mass)
      }
      masses <- sapply(peptide_vector, residueMass)
      pm <- sum(masses)
      p1 <- round(pm + H * 2 + O + proton, digits = 5)
      if (fragments == "by") {
        b1 <- vector(mode = "numeric", length = 0)
        bi <- vector(mode = "integer", length = 0)
        y1 <- vector(mode = "numeric", length = 0)
        yi <- vector(mode = "integer", length = 0)
        for (i in 1:(peptide_length - 1)) {
          mass <- sum(masses[1:i])
          b1[i] <- round(mass + proton, digits = 5)
          bi[i] <- i
        }
        for (j in 2:peptide_length) {
          mass <- sum(masses[j:peptide_length])
          y1[j - 1] <- round(mass + H * 2 + O + proton, digits = 5)
          yi[j - 1] <- peptide_length - j + 1
        }
        ms1z1 <- rep(p1, times = (length(bi) + length(yi)))
        b1.type <- paste("b", bi, sep = "")
        y1.type <- paste("y", yi, sep = "")
        ms2type <- c(b1.type, y1.type)
        MH <- c(b1,y1)
      }
      if (fragments == "cz") {
        c1 <- vector(mode = "numeric", length = 0)
        ci <- vector(mode = "integer", length = 0)
        z1 <- vector(mode = "numeric", length = 0)
        zi <- vector(mode = "integer", length = 0)
        for (i in 1:(peptide_length - 1)) {
          mass <- sum(masses[1:i])
          c1[i] <- round(mass + 3 * H + N + proton, digits = 5)
          ci[i] <- i
        }
        for (j in 2:peptide_length) {
          mass <- sum(masses[j:peptide_length])
          z1[j - 1] <- round(mass + O - N + proton, digits = 5)
          zi[j - 1] <- peptide_length - j + 1
        }
        ms1z1 <- rep(p1, times = (length(ci) +length(zi)))
        c1.type <- paste("c", ci, sep = "")
        z1.type <- paste("z", zi, sep = "")
        ms2type <- c(c1.type, z1.type)
        MH <- c(c1, z1)
      }
      results_list[[sequence_number]] <- data.frame(ms1z1,Ion = ms2type,ms2type, MH)%>% tidyr::separate(ms2type, c("Ion_type", "Ion_num"), sep = 1)
    }
    return(as.data.frame(do.call("rbind", results_list)))
  }
  #ref of WT
  WT_ref <- dplyr::select(WT_ref, c(Fragments,	Ion.type, Ion.num))
  #PredictDiag
  WT <- WT$`sp|P68871|HBB_HUMAN`
  f <- list()
  for (i in 1:length(names(Hbvariants))){
    MT <- matrix(Hbvariants[[i]], byrow = TRUE)
    z <- MT==WT
    f[[i]]<- data.frame(variant=names(Hbvariants)[i], mut.site_N = which(z==FALSE)-1,mut.site_C=148-which(z==FALSE),
                         WT.AA= WT[which(z==FALSE)],MT.AA= MT[which(z==FALSE)])
  }
  f1 <- do.call(rbind,f)
  f1$Mutation <- paste(f1$WT.AA, f1$mut.site_N,f1$MT.AA)
  IDs <- f1$mut.site_N
  mutMass <- function(residue)
  { if (residue == "A")
    mass = 71.03711
  if (residue == "R")
    mass = 156.10111
  if (residue == "N")
    mass = 114.04293
  if (residue == "D")
    mass = 115.02694
  if (residue == "E")
    mass = 129.04259
  if (residue == "Q")
    mass = 128.05858
  if (residue == "G")
    mass = 57.02146
  if (residue == "H")
    mass = 137.05891
  if (residue == "I")
    mass = 113.08406
  if (residue == "L")
    mass = 113.08406
  if (residue == "K")
    mass = 128.09496
  if (residue == "M")
    mass = 131.04049
  if (residue == "F")
    mass = 147.06841
  if (residue == "P")
    mass = 97.05276
  if (residue == "S")
    mass = 87.03203
  if (residue == "T")
    mass = 101.04768
  if (residue == "W")
    mass = 186.07931
  if (residue == "Y")
    mass = 163.06333
  if (residue == "V")
    mass = 99.06841
  if (residue == "C")
    mass = 103.00919
  return(mass)
  }
  f1$mutMass <- round(sapply(f1$MT.AA,mutMass)-sapply(f1$WT.AA,mutMass), 5)

  df <- dplyr::left_join(f1, diag_ref, by = "mut.site_N")
  df2 <- df[, c(1,2,3,6,7,9,10,11,12)]

  #2
  #bc ions
  HbA.bc1 <- subset(WT_ref, Ion.type=="b"|Ion.type=="c")
  HbA.bc <- data.table::as.data.table(HbA.bc1)
  #yz ions
  HbA.yz1 <- subset(WT_ref, Ion.type=="y"|Ion.type=="z")
  HbA.yz <-data.table::as.data.table(HbA.yz1)
  #bc ions
  m <- length(df2$diag.N.term.start)
  f <- list()
  for (i in 1:m){
    s <- as.numeric(df2$diag.N.term.start[i])
    e <- as.numeric(df2$diag.N.term.end[i])
    f[[i]] <- cbind(variant=df2$variant[i],HbA.bc[between(Ion.num,s, e)], Mutation = df2$Mutation[i], MutMass.Da= df2$mutMass[i])
  }
  f2 <- do.call(rbind,f)
  #yz ions
  n <- length(df2$diag.C.term.start)
  f3 <- list()
  for (j in 1:n){
    s <- as.numeric(df2$diag.C.term.start[j])
    e <- as.numeric(df2$diag.C.term.end[j])
    f3[[j]] <- cbind(variant=df2$variant[j],HbA.yz[between(Ion.num,s, e)],Mutation = df2$Mutation[j], MutMass.Da= df2$mutMass[j])
  }
  f4 <- do.call(rbind,f3)

  f5 <- rbind(f2, f4) %>% dplyr::filter(!is.na(Fragments))
  f6 <- f5[with(f5, order(variant))]
  #get monomass
  ID2 <- names(Hbvariants)
  out2 <- list()
  for (j in 1:length(ID2)){
    Hbseq1 <- seqinr::getSequence(Hbvariants, as.string = TRUE)
    Hbseq2 <- sub("M","",Hbseq1[[j]])
    input <- monomz(Hbseq2, fragments = "by")

    input2 <- monomz(Hbseq2, fragments = "cz")

    out1 <- rbind(input,input2)
    out1$varints <- ID2[j]
    out2[[j]] <- out1
  }
  all.frags <- do.call(rbind,out2)

  ID3 <- unique(f6$variant)
  out <- list()
  for (i in 1:length(ID3)){
    sub <- subset(f6, variant==ID3[i])
    sub %>% dplyr::mutate_if(is.numeric, as.character) -> sub
    sub2 <- subset(all.frags, varints==ID3[i])
    sub2 %>% dplyr::mutate_if(is.numeric, as.character) -> sub2
    out[[i]] <- dplyr::inner_join(sub2, sub, by = c("Ion_type"="Ion.type", "Ion_num"="Ion.num"))
  }
  flist <- do.call(rbind,out)
  relist2 <- relist %>% tidyr::separate(variant,c("x","Variant"), sep = 3)
  names(relist2)[5] <- "Ref_Mass"
  relist2$Name <- paste(relist2$Ion,relist2$Variant, sep = "_" )
  relist3 <- relist2[, c(12,2,3,4,5,8,10,11)]
  return(relist3)
}
