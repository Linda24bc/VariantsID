Variants.Identifier <- function(ref, exp, ppm_error_start=-2, ppm_error_end=5){
  exp <- mutate(exp,
                Exp_rel_abundance = Exp_Intensity/max(exp$Exp_Intensity)*100)
  name.ref <- as.character(ref$Name)
  exp$Name <- numeric(length(exp$Exp_Mass))
  for (i in 1:length(exp$Exp_Mass)){
    for (j in 1:length(ref$Ref_Mass)){
      ppm_error = (exp$Exp_Mass[i] - ref$Ref_Mass[j])/ref$Ref_Mass[j] * (10^6)
      if (ppm_error <= ppm_error_end & ppm_error > ppm_error_start)
        exp$Name[i]<- name.ref[j]
    }
  }
  join <- full_join(exp, ref)
  join$ppm_error=(join$Exp_Mass - join$Ref_Mass)/join$Ref_Mass * (10^6)
  join1 <- filter(join, Exp_rel_abundance > 2 & !is.na(Ref_Mass) & !is.na(Exp_Mass))
  sort <- join1[with(join1, order(Variant, Ion.type, Exp_Mass, Ref_Mass)),]
  re <- sort[,c(4,1,3,5,6,7,8,9,10,11,12)]
}
