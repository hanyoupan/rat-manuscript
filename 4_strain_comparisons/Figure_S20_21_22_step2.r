
chrList = c(1:20, "X","Y")

##############################################################
## Merge chromosomes together.
##############################################################

first = TRUE
for (chr in chrList){
  print(paste("adding", chr))
  rolled_overview = read.csv(paste("data/deepvariant_chr",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HIGH_IMPACT_only.gvcf.rolled.gt",sep = ""), sep = "\t", header = T)
  high_impact = read.csv(paste("data/deepvariant_chr",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HIGH_IMPACT_only.gvcf.high.gt",sep = ""), sep = "\t", header = T)
  
  if (first == TRUE){
      all_high = high_impact
      allRolled = rolled_overview
      first = FALSE
  }else{
      all_high = rbind(all_high, high_impact)
      allRolled = rbind(allRolled, rolled_overview)
  }
}


selection = (all_high$ACI.N_JL.ILM == "1/1" & all_high$impact == "HIGH")
overview = all_high[selection,]

all_high_binary = ifelse(all_high[,9:171] == "1/1", 1,0)

################################################################################
# Answering Questions
################################################################################

length(unique(all_high$ID))
# There are 5879 high impact SNPs
length(unique(all_high$nearby_gene))
# Nearby 4431 unique genes
sort(colSums(all_high_binary))
# WKY.NCrl_MD.ILM has 1305 High impact SNPs.

write.csv(sort(table(all_high$region_type), decreasing = T),"results/01_highImpact_regionType.csv")

length(unique(all_high$ID))

chromosomalMax = tapply(allRolled$end, allRolled$chromosome, max)

strainName = "ACI.N_AP.ILM"
localColor = c("black","darkgrey")

write.csv(all_high, "data/high_impact_snps.csv")

