#install.packages("tidyverse")
#install.packages("Polychrome")

library(forcats)
library(scales)
library(tidyverse)
library(magrittr)
library('rlist')
library("gplots")

numberToPlot = 20
yAxisSize = 0.9
xAxisSize = 0.5

#remove.packages('gplots'); 
#library('devtools'); 
#install_github("ChristophH/gplots")
library("gplots")

library("Polychrome")

myColors = c('#e6194b', '#4363d8' , '#ffe119','#3cb44b', 
             '#f58231', '#911eb4', '#46f0f0', '#f032e6',
             '#bcf60c', '#fabebe', '#008080', '#e6beff',
             '#9a6324', 'gold', '#800000', '#aaffc3',
             '#808000', 'black', '#000075', '#808080')

localColors = c("brown","#6CA2EA","#442288", "#B5D33D","#FED23F")

ratColors = c("brown","lightblue","#651B20","orange","pink","ivory","#5D594F","#C8AC97")


library("Polychrome")
data(palette36)
newColors = palette36


strainOrder = read.csv(file = "name2_ibs_order.txt", header = F)
strainOrder$V1 = gsub("\\*","", strainOrder$V1)

ibsOrder = read.csv(file = "ibs_order.txt", header = F)


################################################################################
# 1. Set chromosome sizes
################################################################################

chrList = c(1:20,"X","Y")

chromosome_sizes = c(282763074, # 1
                     266435125, # 2
                     175699992, # 3
                     184226339, # 4
                     173707219, # 5
                     147991367, # 6
                     145729302, # 7
                     133307652, # 8
                     122095297, # 9
                     112626471, # 10
                     90463843,  # 11
                     52716770,  # 12
                     114033958, # 13
                     115493446, # 14
                     111246239, # 15
                     90668790,  # 16
                     90843779,  # 17
                     88201929,  # 18
                     62275575,  # 19
                     56205956,  # 20
                     159970021, # X
                     3310458)   # Y

chromosome_sizes = chromosome_sizes 

chromosomal_addition = 0 
chromosome_backToBack = c(0)

for(i in chromosome_sizes){
  chromosomal_addition = chromosomal_addition + i + 20000000
  print(chromosomal_addition)
  chromosome_backToBack = c(chromosome_backToBack, (chromosomal_addition))
}
chromosome_backToBack = chromosome_backToBack 


################################################################################
# 2. Getting associated disease information
################################################################################

# Oke, we have the disease data. Now I need a way to link each gene to the disease data and color code it.
# There are multiple disease triats, perhaps sort these first?

inData = read.delim("data/high_impact_annotated.csv", sep = "\t")
dim(inData)
head(inData)
table(inData$impact)

colnames(inData)[colnames(inData) == "FXLE18_MD..ILM"] = "FXLE18_MD.ILM" 


################################################################
## 3. Split the data and put it back together again
################################################################
dataOrder = read.delim("r163_order.txt", sep = "")
newOrder = gsub("/",".",dataOrder[,1])
newOrder = gsub("-",".",newOrder)
################################################################
lastSample = 172

part1 = 1:9
part2 = 10:lastSample
part3 = 173:ncol(inData)

inData_part1 = inData[,part1]
inData_part2 = inData[,part2]
inData_part3 = inData[,part3]

inData = cbind(inData_part1, inData_part2, inData_part3)


################################################################################
# 4. Simplify disease terms
################################################################################

# Lets see the number of diseases:
dim(table(inData$TERM_NAME))

# Perhaps I can broadly group them somehow
inData$simpleTerms = tolower(inData$TERM_NAME)


# 1 cancer
inData$simpleTerms = gsub(".*cancer.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*carcinoma.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*tumor.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*lynch syndrome.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*malignant.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*melanoma.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*carcino.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*xanthoastrocytoma.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*neoplasm.*", "cancer/tumor", inData$simpleTerms)
inData$simpleTerms = gsub(".*neoplastic.*", "cancer/tumor", inData$simpleTerms)


inData$simpleTerms = gsub(".*epilepsy.*", "epilepsy", inData$simpleTerms)

# 2 blood
inData$simpleTerms = gsub(".*blood.*", "blood", inData$simpleTerms)
inData$simpleTerms = gsub(".*hypertens.*", "hypertension", inData$simpleTerms)

# 3 Heart
inData$simpleTerms = gsub(".*Heart.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*cardiac.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*cardio.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*cardial.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*artery.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*aortic.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*atrial.*", "cardiovascular", inData$simpleTerms)
inData$simpleTerms = gsub(".*ventricular.*", "cardiovascular", inData$simpleTerms)

inData$simpleTerms = gsub(".*brugada syndrome.*", "cardiovascular", inData$simpleTerms)

# 4 bladder
inData$simpleTerms = gsub(".*bladder.*", "bladder", inData$simpleTerms)

# 5 brain
inData$simpleTerms = gsub(".*brain.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*neuro.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*schizophrenia.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*body dysmorphic.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*depression.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*adhd.*", "psychiatric", inData$simpleTerms)
inData$simpleTerms = gsub(".*addiction.*", "psychiatric", inData$simpleTerms)

inData$simpleTerms = gsub(".*epileptic.*", "epilepsy", inData$simpleTerms)
inData$simpleTerms = gsub(".*epilepsy.*", "epilepsy", inData$simpleTerms)


# 6 lung
inData$simpleTerms = gsub(".*lung.*", "lung associated", inData$simpleTerms)
inData$simpleTerms = gsub(".*asthma.*", "asthma", inData$simpleTerms)
inData$simpleTerms = gsub(".*COPD.*", "COPD", inData$simpleTerms)
inData$simpleTerms = gsub(".*pulmonary.*", "lung associated", inData$simpleTerms)
inData$simpleTerms = gsub(".*respiratory.*", "lung associated", inData$simpleTerms)
#inData$simpleTerms = gsub(".*covid-19.*", "lung", inData$simpleTerms)
inData$simpleTerms = gsub(".*pneumo.*", "lung associated", inData$simpleTerms)

# 7 liver
inData$simpleTerms = gsub(".*liver.*", "liver", inData$simpleTerms)

# 8 kidney
inData$simpleTerms = gsub(".*kidney.*", "kidney", inData$simpleTerms)

# 9 immune system
inData$simpleTerms = gsub(".*immuno.*", "immune", inData$simpleTerms)
inData$simpleTerms = gsub(".*immune.*", "immune", inData$simpleTerms)

# 10 arthritis
inData$simpleTerms = gsub(".*arthritis.*", "arthritis", inData$simpleTerms)

# 11 bone
inData$simpleTerms = gsub(".*bone.*", "bone", inData$simpleTerms)

# 12 skin
inData$simpleTerms = gsub(".*skin.*", "skin", inData$simpleTerms)

# 13 obesity
inData$simpleTerms = gsub(".*obesity.*", "obesity", inData$simpleTerms)

# 14
inData$simpleTerms = gsub(".*insulin.*", "diabetes", inData$simpleTerms)
inData$simpleTerms = gsub(".*diabetic.*", "diabetes", inData$simpleTerms)

# 15
inData$simpleTerms = gsub(".*intellectual disability.*", "intellectual disability", inData$simpleTerms)
inData$simpleTerms = gsub(".*mental retardation.*", "intellectual disability", inData$simpleTerms)

# 16
inData$simpleTerms = gsub(".*deafness.*", "deafness", inData$simpleTerms)

# 17
inData$simpleTerms = gsub(".*muscle.*", "muscle", inData$simpleTerms)
inData$simpleTerms = gsub(".*muscular.*", "muscle", inData$simpleTerms)
inData$simpleTerms = gsub(".*myofibrillar.*", "muscle", inData$simpleTerms)
inData$simpleTerms = gsub(".*myopathy.*", "muscle", inData$simpleTerms)

inData$simpleTerms = gsub(".*retinitis.*", "retinitis", inData$simpleTerms)


# 18
inData$simpleTerms = gsub(".*diabetes.*", "diabetes", inData$simpleTerms)

# 19
inData$simpleTerms = gsub(".*autism.*", "autism spectrum disorder", inData$simpleTerms)
inData$simpleTerms = gsub(".*autistic.*", "autism spectrum disorder", inData$simpleTerms)


table(inData$simpleTerms)
dim(table(inData$simpleTerms))

head(sort(table(inData$simpleTerms),decreasing = T))

stored_inData = inData


################################################################################
# 5. How many annotated high impact snps?
################################################################################

#uniqueData = inData[!duplicated(inData$ID),]
uniqueData = inData[!duplicated(inData$nearby_gene),]

# Number of High impact SNPs with annotation
length(unique(inData$ID))

# Number of unique genes affected
length(unique(uniqueData$nearby_gene))

# Number of unique disease terms
length(unique(inData$TERM_NAME))

# Strain with most high impact SNPs
sort(colSums(uniqueData[,part2] == "1/1")) 

# Strain with most high impact SNPs
sort(colSums(uniqueData[,part2] == "1/1"))

# Number of disease terms
length(unique(inData$TERM_NAME))

head(sort(table(inData$simpleTerms),decreasing = T), n = 20)


# Get the unique genes first. Then use these names to get disease terms
unique(inData$ID)
# There are 2601 high impact SNPs


#########################################
## 6. What are the most common region types
#########################################
uniqueData = inData[!duplicated(inData$ID),]
regionTypes = table(uniqueData$region_type)
regionTypes = regionTypes[order(regionTypes,decreasing = T)]
write.csv(regionTypes, "03_regionType_HighImpact.csv")


#########################################
## 7. What are the most common region types
#########################################

print(length(unique(inData$simpleTerms)))
first = TRUE
for (simple in unique(inData$simpleTerms)){
  print(simple)
  if(first == TRUE){
    first = FALSE
    fullNames = unique(inData$TERM_NAME[inData$simple == simple])
    mergedListOverview = cbind(fullNames, rep(simple, length(fullNames)))
    colnames(mergedListOverview) = c("Simplfied disease name","Full disease name")
    
  }else{
    fullNames = unique(inData$TERM_NAME[inData$simple == simple])
    simpleListOverview = cbind(fullNames, rep(simple, length(fullNames)))
    mergedListOverview = rbind(mergedListOverview, simpleListOverview)
  }
}

head(mergedListOverview)
write.csv(mergedListOverview, "results/01_simplified_disease_name_index.csv")

inData$limited_groups = paste(inData$ID,inData$simpleTerms,sep = "_")
inData = inData[!duplicated(inData$limited_groups),]


#############################################################################################################
# 8. Simplify sample names
#############################################################################################################

# Create plottable groups
simplifiedNames = colnames(inData)[part2]
simplifiedNames =  gsub("*ACI.*", "ACI", simplifiedNames)

simplifiedNames =  gsub("BBDP.*", "BBDP", simplifiedNames)
simplifiedNames =  gsub(".*BN\\.Lx.*" , "BN_LX",  simplifiedNames)
simplifiedNames =  gsub(".*BN\\..*" , "BN",  simplifiedNames)
simplifiedNames =  gsub("BUF.*", "BUF", simplifiedNames)
#simplifiedNames =  gsub("BXH.*", "BXH", simplifiedNames)

simplifiedNames =  gsub("DA.*", "DA", simplifiedNames)

simplifiedNames =  gsub("FHH.*", "FHH", simplifiedNames)
simplifiedNames =  gsub("FHL.*", "FHL", simplifiedNames)
simplifiedNames =  gsub("F344.*", "F344", simplifiedNames)

simplifiedNames =  gsub("GK.*", "GK", simplifiedNames)

simplifiedNames =  gsub("HSRA.*", "HSRA", simplifiedNames)
#simplifiedNames =  gsub("HXB.*", "HXB", simplifiedNames)

simplifiedNames =  gsub("LL.*", "LL", simplifiedNames)
simplifiedNames =  gsub("LH.*", "LH", simplifiedNames)
simplifiedNames =  gsub("LN.*", "LN", simplifiedNames)

simplifiedNames =  gsub("LEW.*", "LEW", simplifiedNames)
#simplifiedNames =  gsub("LEXF.*", "LEXF", simplifiedNames)
#simplifiedNames =  gsub(".*LXF10A.*", "LEXF", simplifiedNames)
#simplifiedNames =  gsub("FXLE.*", "FXLE", simplifiedNames)

simplifiedNames =  gsub("M520.*", "M520", simplifiedNames)
simplifiedNames =  gsub("MR.*", "MR", simplifiedNames)

simplifiedNames =  gsub("MHF.*", "MHF", simplifiedNames)
simplifiedNames =  gsub("MHS.*", "MHS", simplifiedNames)
simplifiedNames =  gsub("MNS.*", "MNS", simplifiedNames)
simplifiedNames =  gsub("MWF.*", "MWF", simplifiedNames)

simplifiedNames =  gsub("PVG.*", "PVG", simplifiedNames)

simplifiedNames =  gsub(".*SBH.*", "SBH", simplifiedNames)
simplifiedNames =  gsub(".*SBN.*", "SBN", simplifiedNames)
simplifiedNames =  gsub("SHR.*", "SHR", simplifiedNames)
simplifiedNames =  gsub("SS.*", "SS", simplifiedNames)
simplifiedNames =  gsub("SR.*", "SR", simplifiedNames)

simplifiedNames =  gsub("WAG.*", "WAG", simplifiedNames)
simplifiedNames =  gsub("WKY.*", "WKY", simplifiedNames)
simplifiedNames =  gsub("WLI.*", "WLI", simplifiedNames)
simplifiedNames =  gsub("WMI.*", "WMI", simplifiedNames)
simplifiedNames =  gsub("WN.*", "WN", simplifiedNames)
simplifiedNames =  gsub("LE.Stm.*", "LE", simplifiedNames)
simplifiedNames =  gsub("ACI.*", "ACI", simplifiedNames)

sampleFactors = as.factor(simplifiedNames)
sampleColors = c(newColors[sampleFactors])

# Check
table(simplifiedNames == colnames(inData)[10:lastSample])

length(table(simplifiedNames))

# Lets write the tool to plot horizontally.
chromosome_sizes
chromosome_backToBack


inData_backup = inData

for(i in part2){
  
  inData[,i] = ifelse(inData[,i] == "1/1",1,0)
}

sampleNames = unique(simplifiedNames)
numberOfSamples = length(sampleNames)

interestingNames = c()
for (simple in simplifiedNames){
  print (simple)
  simpleData = inData[part2][,simplifiedNames == simple]
  
  if (sum(simplifiedNames == simple)>1){
    inData[,simple] = rowSums(simpleData == 1)
  }else{
    inData[,simple] = simpleData
  }
}

################################################################################
## 9. Per larger group
################################################################################

# Select only the simplified strain names
head(inData[,191:ncol(inData)])
StrainData = ifelse(inData[191:ncol(inData)] > 0, 1,0)

allData = inData[]

# Only get the top 20 diseases
diseaseInfo = inData$simpleTerms
LabelsOfInterest = head(sort(table(diseaseInfo),decreasing = TRUE), n = 50)
filter = diseaseInfo %in% names(LabelsOfInterest)

# Apply the filter
plotData = StrainData[filter,]
diseaseInfo = diseaseInfo[filter]

write.csv(LabelsOfInterest, file = "results/top_disease_terms.csv")

################################################################################
## 10. In the phylogenic order
################################################################################


## Load a table with the High impact SNPs
sortedData = inData
sortedData[,part2] = ifelse(sortedData[,part2] > 0, 1, 0)

 # Create filter for top 20 diseases
diseaseInfo = inData_backup$simpleTerms
LabelsOfInterest = head(sort(table(diseaseInfo),decreasing = TRUE), n = numberToPlot)

filter = diseaseInfo %in% names(LabelsOfInterest) & (rowSums(sortedData[,part2]) > 1)

plotInfo = diseaseInfo[filter]
plotData = sortedData[filter,][,part2]

plotOrder = order(plotInfo,-rowSums(plotData))

par(mar = c(20.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right)

simpleNames = gsub("_.*$","",colnames(plotData))

# Create labels for the rat strains
shortNames = substr(simpleNames,1,3)
shortNames = sub("BXH","HXB",shortNames)
shortNames = sub("LEXF","FXLE",shortNames)
shortColors = rainbow(36)
table(shortNames)
length(table(shortNames))


colorsOfInterest = rainbow(n=10)[as.factor(plotInfo)]

plotColors = c("#9FE3FE", "#9CB5FE",
               "#DF9BFD", "#FD9DCB",
               "#FD9DCB", "#FECF9F",
               "#FEE9A7", "#F7FCA8",
               "#D3FFA5","forestgreen","black")[as.factor(plotInfo)]

dim(plotData[plotOrder,])


#boxData = boxplot(plotData ~ plotInfo, las = 2)
#boxData


## Original
heatmap.2(as.matrix(plotData[plotOrder,]),
          dendrogram = "none",
          labCol = simpleNames,
          labRow = plotInfo[plotOrder],
          cexCol = yAxisSize,
          cexRow = 1.0,
          scale = "none",
          Colv = FALSE,
          Rowv = FALSE,
          col = c("white","red"),
          trace = "none",
          margins = c(10,12),
          RowSideColors = plotColors[plotOrder],
          ColSideColors = shortColors[as.factor(shortNames)])

aggregated_Data = aggregate(.~plotInfo, data = as.data.frame(plotData), sum)

rownames(aggregated_Data) = aggregated_Data[,1]

head(aggregated_Data)
aggregated_Data = aggregated_Data[,-1]


library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

## Aggregated
heatmap.2(as.matrix(aggregated_Data),
          dendrogram = "none",
          #labCol = simpleNames,
          #labRow = plotInfo[plotOrder],
          cexCol = 1,
          #cexRow = 1.0,
          scale = "row",
          Colv = FALSE,
          Rowv = FALSE,
          col = Colors,
          trace = "none",
          margins = c(10,10),
          
          #RowSideColors = plotColors[plotOrder],
          ColSideColors = shortColors[as.factor(shortNames)])

## Lets group the samples in a different manner:

groupNames = colnames(aggregated_Data)
shortNames = substr(groupNames,1,3)

founders = c("ACI.N_JL.ILM","ACI.N_AP.ILM",
             "BN.N_JL.ILM","BN.N_AP.ILM",
             "BUF.N_JL.ILM","BUF.N_AP.ILM",
             "F344.N_JL.ILM","F344.N_AP.ILM",
             "M520.N_JL.ILM","M520.N_AP.ILM",
             "MR.N_JL.ILM","MR.N_AP.ILM",
             "WN.N_JL.ILM","WKY.N_AP.ILM",
             "WKY.N_JL.ILM","WN.N_AP.ILM")

BNLX = c("BN.Lx.Cub_HCJL.CRM" ,"BN.Lx.Cub_MD.ILM")

groupNames = ifelse(groupNames %in% founders, "HS Founders", "Classic Inbreds")
groupNames = ifelse(shortNames == "BN.", "BN (Reference strain)", groupNames)
groupNames = ifelse(shortNames == "HXB" | shortNames == "BXH", "HXB/BXH", groupNames)
groupNames = ifelse(colnames(aggregated_Data) %in% BNLX, "HXB/BXH", groupNames)

groupNames = ifelse(shortNames == "LEX" | shortNames == "FXL", "FXLE/LEXF", groupNames)
groupNames = ifelse(shortNames == "LXF", "FXLE/LEXF", groupNames)

test1 = cbind(groupNames, colnames(aggregated_Data))

groupNames = factor(groupNames)

table(groupNames)

# founders, FXLE, HXB,Reference, Selectively bred
displayNames = sub("_.*", "",colnames(aggregated_Data))





## Do this once



################################################################################
## 11. For publication
################################################################################
intersect(colnames(aggregated_Data),newOrder)
setdiff(colnames(aggregated_Data),newOrder)
setdiff(newOrder,colnames(aggregated_Data))
aggregated_Data = aggregated_Data[,newOrder]



# How many unique disease terms?
head(plotInfo)
head(plotData)

## Figure with the top 20 disease terms:

pdf("figure_diseasePerSample.pdf", width=12, height=6)
## Aggregated
heatmap.2(as.matrix(aggregated_Data),
          dendrogram = "none",
          Colv = FALSE,
          Rowv = FALSE,
          labCol = displayNames,
          #labRow = plotInfo[plotOrder],
          offsetCol = 0.5,
          cexCol = xAxisSize,
          cexRow = yAxisSize,
          scale = "row",
          col = Colors,
          trace = "none",
          margins = c(10,12),
          density.info="none",
          keysize = 2,
          key.xlab = "Relative number of SNPs in disease category\n(Z-score)",
          #RowSideColors = plotColors[plotOrder],
          ColSideColors = localColors[groupNames])
legend(x = "topright", 
       legend = names(table(groupNames)),
       col = localColors, 
       pch = 15, cex = 1, pt.cex = 2,
       ncol = 1,
       inset = c(-0.0, -0.0),
       xpd = T, bty = "n",
       horiz = F
)


dev.off()


aggregated_Data = aggregate(.~plotInfo, data = as.data.frame(plotData), sum)

head(inData[,191:ncol(inData)])

collapsed = inData[,191:ncol(inData)] 

my_info = inData$simpleTerms


# Should be done per column # Edit here for Fractions.
collapsed_step = sapply(collapsed, max, na.rm = TRUE)
collapsed_fraction = collapsed
for(i in 1:ncol(collapsed)){
  print(i)
  colnames(collapsed_fraction)[i]
  collapsed_fraction[,i] = collapsed_fraction[,i] / collapsed_step[i]
}

collapsed_binary = ifelse(collapsed >= 1,1,0)

#collapsed_binary = collapsed / 



##########################################################################################################################
# B. Final figure Fractional version. (A. was deleted).
##########################################################################################################################

short_names = substr(colnames(inData),1,3)
HXB_list = short_names %in% c("HXB", "BXH", "FXL","LEX","LXF")
allTheHXBs = inData[,HXB_list]

collapsed_binary2 = cbind(collapsed_binary,allTheHXBs)
collapsed_fraction2 = cbind(collapsed_fraction,allTheHXBs)

aggregated_Data_collapsed = aggregate(.~my_info, data = as.data.frame(collapsed_fraction2), sum) # Collapsed version


################################################################################

#               founders, FXLE,   HXB,        Reference, Selectively bred

plot(c(1,2,3,4,5), col = localColors, pch =15)

test = aggregated_Data_collapsed[order(-rowSums(aggregated_Data_collapsed[,2:ncol(aggregated_Data_collapsed)])),]
rownames(test) = test[,1]
test = test[,-1]

testPart1 = test[,colnames(collapsed_fraction)]
testPart2 = test[,colnames(allTheHXBs)]

testPart1 = testPart1[,order(-colSums(testPart1))]
testPart2 = testPart2[,order(-colSums(testPart2))]
colnames(testPart2) = gsub("_.*", "", colnames(testPart2))

test = cbind(testPart1,testPart2)
colnames(test) = gsub("\\.1.*", "",colnames(test))
 

founders = c("ACI","BUF","F34","M52","MR","WN","WKY")
shortNames = substr(colnames(test),1,3)
shortNames = ifelse(shortNames %in% founders, "founder",shortNames)

newColors = ifelse(shortNames == "FXL" | shortNames == "LXF" | shortNames == "LEX", "#442288", "#6CA2EA") # FXL, Others
newColors = ifelse(shortNames == "HXB" | shortNames == "BXH" , "#FED23F", newColors) # HXB
newColors = ifelse(shortNames == "BN_" | shortNames == "BN"  , "brown", newColors) # BN
newColors = ifelse(shortNames == "founder" , "#B5D33D", newColors) # Founders

localColors = c("brown","#6CA2EA","#B5D33D", "#442288" ,"#FED23F")
localNames = c("BN (Reference strain)","Classic Inbreds","HS Founders", "FXLE/LEXF","HXB/BXH" )  

################################################################################
## Top 5 disease ontology per gene and per snp.
################################################################################

testFilter_terms_gene = paste(inData$simpleTerms,inData$nearby_gene)
testFilter_terms_ID = paste(inData$simpleTerms,inData$ID)

inData_termGene = inData[!duplicated(testFilter_terms_gene),]
inData_termID = inData[!duplicated(testFilter_terms_ID),]

# 1. How many Genes with annotation
length(unique(inData$nearby_gene)) # 2079 unique genes with annotation

# 2.  How many unique annotations
length(unique(inData$TERM_NAME)) # 3268 unqiue disease terms

# 3. What are the top 5, or top 10 disease terms
head(sort(table(inData_termGene$simpleTerms), decreasing = T), n = 10)

head(sort(table(inData_termID$simpleTerms), decreasing = T), n = numberToPlot)
plotThese = head(sort(table(inData_termID$simpleTerms), decreasing = T), n = numberToPlot)

################################################################################
## Version 1. Everything dendrogrammed.
################################################################################
pdf("02_disease_overview_fractioned.pdf", width=10, height=5)

## Aggregated
plotMe = test[1:numberToPlot,][order(rownames(test[1:numberToPlot,])),]
plotMe = test[names(plotThese),]
plotLabels1 = paste(names(plotThese)," (",plotThese, ")", sep = "")
plotLabels2 = paste(names(plotThese), "")


################################################################################
## Version 2. Manual sorted
################################################################################
## Aggregated

# Ignore this lazy fix
colnames(plotMe) = ifelse(colnames(plotMe) =="HSR", "HSRA",colnames(plotMe))

plotMe = round(plotMe, 0)

#################
## Version 1
#################

heatmap.2(as.matrix(plotMe),
          dendrogram = "none",
          Colv = FALSE,
          Rowv = FALSE,
          #labCol = displayNames,
          labRow = plotLabels1,
          offsetCol = 0.5,
          cexCol = 0.6,
          cexRow = yAxisSize,
          scale = "row",
          col = Colors,
          trace = "none",
          margins = c(10,12),
          density.info="none",
          keysize = 1,#,
          key.xlab = "Relative number of SNPs in disease category\n(Z-score)",
          #RowSideColors = plotColors[plotOrder],
          ColSideColors = newColors)

groupNames

legend(x = "topright", 
       legend = localNames,
       col = localColors, 
       pch = 15, cex = 1, pt.cex = 2,
       ncol = 1,
       inset = c(-0.0, -0.0),
       xpd = T, bty = "n",
       horiz = F
)

#################
## Version 2
#################

heatmap.2(as.matrix(plotMe),
          dendrogram = "none",
          Colv = FALSE,
          Rowv = FALSE,
          #labCol = displayNames,
          labRow = plotLabels2,
          offsetCol = 0.5,
          cexCol = 0.6,
          cexRow = yAxisSize,
          scale = "row",
          col = Colors,
          trace = "none",
          margins = c(10,12),
          cellnote = plotMe,
          notecex=0.3,
          notecol = "black",
          density.info="none",
          keysize = 2,#,
          key.xlab = "Relative number of SNPs in disease category\n(Z-score)",
          #RowSideColors = plotColors[plotOrder],
          ColSideColors = newColors)

legend(x = "topright", 
       legend = localNames,
       col = localColors, 
       pch = 15, cex = 1, pt.cex = 2,
       ncol = 1,
       inset = c(-0.0, -0.0),
       xpd = T, bty = "n",
       horiz = F
)

dev.off()

################################################################################
## HXB only, LEXF only, other strain only
################################################################################

makeHeatmap = function(plotMe, totalSNP, xlabSize, fileName,bigPlot){
  testLabel = ""
  selection = 1:numberToPlot
  outDisease = paste(rownames(plotMe[selection,])," (",totalSNP,")",sep = "")
  outStrain = sub("\\.1", "",colnames(plotMe))
  
  plotMe = round(plotMe,0)
  
  pdf(paste0("results/03_",fileName, ".pdf"),width = (ncol(plotMe)+5)/4, height = nrow(plotMe)/4)
  
  if (bigPlot == TRUE){
    extraCol = 0
    extraRow = 0
  }else{
    extraCol = 33 - ncol(plotMe)
    extraRow = 33- nrow(plotMe)
  } 
  
  #squaring the plot
  emptyData = data.frame(matrix(NA, nrow = nrow(plotMe), ncol = extraCol))
  plotMe2 = cbind(emptyData,plotMe)
  outStrain2 = c(rep(testLabel,extraCol),outStrain)
  
  #squaring the plot pt2.
  emptyData2 = data.frame(matrix(NA, nrow = extraRow, ncol = ncol(plotMe2)))
  colnames(emptyData2) = colnames(plotMe2)
  plotMe3 = rbind(emptyData2,plotMe2)
  outDisease2 = c(rep(testLabel,extraRow),outDisease)
  
  # For the big plot
  if (bigPlot == TRUE){

    shortNames = substr(colnames(plotMe3),1,3)
    HXBLEXF_list = shortNames %in% c("HXB", "BXH","LXF","LEX","FXL")
    plotPart1 = plotMe3[,!HXBLEXF_list]
    plotPart2 = plotMe3[,HXBLEXF_list]
    plotMe3 = cbind(plotPart1,plotPart2)
    outStrain2 = colnames(plotMe3)
    
    founders = c("ACI","BN","BUF","F34","M52","MR","WN","WKY")
    shortNames = substr(outStrain2,1,3)
    newColors = ifelse(shortNames %in% c("FXL","LXF", "LEX"), "#442288", "#6CA2EA") # FXL, Others
    newColors = ifelse(shortNames %in% c("HXB","BXH") , "#FED23F", newColors) # HXB
    newColors = ifelse(shortNames %in% c(founders) , "#B5D33D", newColors)
    newColors = ifelse(shortNames == "BN", "brown", newColors)
    # BN
  }else{
    newColors = rep("white",ncol(plotMe3))
  }

  heatmap.2(as.matrix(plotMe3[selection,]),
              dendrogram = "none", #row column
              Rowv=F, Colv=F,
              labCol = outStrain2,
              labRow = outDisease2,
              ColSideColors = newColors,
              offsetCol = 0.0,
              offsetRow = -0.0,
              adjCol = c(NA,0.5),
              cexCol = xlabSize,
              cexRow = yAxisSize,
              scale = "row",
              col = Colors,
              trace = "none",
              margins = c(10,12),
              density.info="none",
              symm = TRUE,
              keysize = 2,#,
              key.xlab = "Relative number of SNPs in disease category\n(Z-score)"
              #RowSideColors = plotColors[plotOrder],
              #ColSideColors = newColors)
              )
  
  if(bigPlot == TRUE){
      legend(x = "topright", 
             legend = localNames,
             col = localColors, 
             pch = 15, cex = 1, pt.cex = 2,
             ncol = 1,
             inset = c(-0.0, -0.2),
             xpd = T, bty = "n",
             horiz = F
      )
  }
 
  heatmap.2(as.matrix(plotMe3[selection,]),
            dendrogram = "none", #row column
            Rowv=F, Colv=F,
            labCol = outStrain2,
            labRow = outDisease2,
            offsetCol = 0.0,
            offsetRow = -0.0,
            adjCol = c(NA,0.5),
            cexCol = xlabSize,
            cexRow = 0.6,
            scale = "row",
            col = Colors,
            trace = "none",
            margins = c(10,12),
            density.info="none",
            symm = TRUE,
            cellnote = plotMe3,
            notecex=0.3,
            notecol = "black",
            keysize = 2,#,
            key.xlab = "Relative number of SNPs in disease category\n(Z-score)"
            #RowSideColors = plotColors[plotOrder],
            #ColSideColors = newColors)
  )  
  dev.off()

}


################################################################################
## Get shortnames from my aggregated set
################################################################################

shortNames = substr(colnames(test),1,3)
HXB_filter = shortNames == "HXB" | shortNames == "BXH" | 
             shortNames == "SHR" | shortNames == "BN_"

LEXF_filter = shortNames == "LEX" | shortNames == "FXL" |  shortNames == "LXF" |
              shortNames == "F34" | shortNames == "LE"

other_filter = (!HXB_filter & !LEXF_filter) | 
                shortNames == "F34" |shortNames == "SHR" | shortNames == "BN" # | shortNames == "LE" 

ALL_filter = rep(TRUE,length(shortNames))


################################################################################
## Get shortnames from my main dataset
################################################################################

checkup = inData[,part2]
noDuplicate_filter = !duplicated(paste(inData$ID,inData$simpleTerms))
checkup = checkup[noDuplicate_filter,]

shortNames_inData = substr(colnames(checkup),1,3)

ALL_filter2 = rep(TRUE, length(shortNames_inData))

HXB_filter2 = shortNames_inData == "HXB" | shortNames_inData == "BXH" | 
  substr(colnames(checkup),1,7) == "SHR.Ola" | substr(colnames(checkup),1,5) == "BN.Lx"

LEXF_filter2 = shortNames_inData == "LEX" | shortNames_inData == "FXL" |  shortNames_inData == "LXF" |
  shortNames_inData == "F34" | shortNames_inData == "LE"

other_filter2 = (!HXB_filter2 & !LEXF_filter2) | 
  shortNames_inData == "F34" | shortNames_inData == "SHR" | shortNames_inData == "BN" # |  shortNames_inData == "LE"


## ALL total SNPs in Panel
localTest_ALL = inData[rowSums(checkup[ALL_filter2]) >= 1,]
ALL_selection = head(sort(table(localTest_ALL$simpleTerms), decreasing = TRUE), n = 30)

## HXB total SNPs in Panel
localTest_HXB = inData[rowSums(checkup[HXB_filter2]) >= 1,]
HXB_selection = head(sort(table(localTest_HXB$simpleTerms), decreasing = TRUE), n = 30)

## LEXF total SNPs in Panel
localTest_LEXF = inData[rowSums(checkup[LEXF_filter2]) >= 1,]
LEXF_selection = head(sort(table(localTest_LEXF$simpleTerms), decreasing = TRUE), n = 30)

## HXB total SNPs in Panel
localTest_other = inData[rowSums(checkup[other_filter2]) >= 1,]
other_selection = head(sort(table(localTest_other$simpleTerms), decreasing = TRUE), n = 30)


# ALL select
ALL_select30 = test[,ALL_filter][names(ALL_selection),]
ALL_short = gsub("\\..*","",colnames(ALL_select30))
ALL_temp = t(ALL_select30)
ALL_select30 = aggregate(.~ALL_short, data = as.data.frame(ALL_temp), mean)
rownames(ALL_select30) = ALL_select30[,1]
ALL_select30 = t(ALL_select30[,-1])

HXB_select30 = test[,HXB_filter][names(HXB_selection),]
HXB_short = gsub("\\..*","",colnames(HXB_select30))
HXB_temp = t(HXB_select30)
HXB_select30 = aggregate(.~HXB_short, data = as.data.frame(HXB_temp), mean)
rownames(HXB_select30) = HXB_select30[,1]
HXB_select30 = t(HXB_select30[,-1])

LEXF_select30 = test[,LEXF_filter][names(LEXF_selection),]
LEXF_short = gsub("\\..*","",colnames(LEXF_select30))
LEXF_temp = t(LEXF_select30)
LEXF_select30 = aggregate(.~LEXF_short, data = as.data.frame(LEXF_temp), mean)
rownames(LEXF_select30) = LEXF_select30[,1]
LEXF_select30 = t(LEXF_select30[,-1])


other_select30 = test[,other_filter][names(other_selection),]
HXB_select30 = HXB_select30[,order(colSums(HXB_select30))]
LEXF_select30 = LEXF_select30[,order(colSums(LEXF_select30))]
other_select30 = other_select30[,order(colSums(other_select30))]
ALL_select30 = ALL_select30[,order(colSums(ALL_select30))]



pdf("results/03_disease_burden_perPanel.pdf")
makeHeatmap(plotMe = HXB_select30, totalSNP = HXB_selection, xlabSize = 0.5, fileName = "HXB",bigPlot = FALSE)
makeHeatmap(plotMe = LEXF_select30, totalSNP = LEXF_selection, xlabSize = 0.5, fileName = "LEXF",bigPlot = FALSE)
makeHeatmap(plotMe = other_select30, totalSNP = other_selection, xlabSize = 0.5, fileName = "other",bigPlot = FALSE)

dev.off()

pdf("results/03_disease_burden_all.pdf", width = 12,height = 5)
makeHeatmap(plotMe = ALL_select30, totalSNP = ALL_selection, xlabSize = 0.5, fileName = "ALL",bigPlot = TRUE)

dev.off()


################################################################################
# Show that SHR is the same in both the HXB and the "other" group
compare1 = HXB_select30[rownames(other_select30),]
compare2 = cbind(compare1$SHR, other_select30$SHR)
rownames(compare2) = rownames(other_select30)
colnames(compare2) = c("HXB SHR", "other SHR")
compare2 
################################################################################
inData

# How many in HXB with annotation?
lookup = inData[rowSums(checkup[HXB_filter2]) >= 1,]
length(unique(lookup$ID))

# How many in FXLE with annotation?
lookup = inData[rowSums(checkup[LEXF_filter2]) >= 1,]
length(unique(lookup$ID))

# Get list, per investigation group, with info, cbound.
shortNames
inData = stored_inData
head(inData)

part1 = 1:9
part2 = 10:172
part3 = 173:ncol(inData)

inData_part1 = inData[,part1]
inData_part2 = inData[,part2]
inData_part3 = inData[,part3]

################################################################################
#Investigations:
################################################################################
stored_inData

shortNames = substr(colnames(inData_part2),1,3)

#Het stock / Founders
founders = c("ACI","BUF","F34","M52","MR.","WN.","WKY")
get2 = inData_part2[,shortNames %in% founders]
selection = rowSums(get2 == "1/1") >= 1
out_Founders = cbind(inData_part1, get2, inData_part3)[selection,]
write.csv(out_Founders, "results/04_highImpact_founders.csv")

# All
get2 = inData_part2
selection = rowSums(get2 == "1/1") >= 1
out_ALL = cbind(inData_part1, get2, inData_part3)[selection,]
write.csv(out_ALL, "results/04_highImpact_All.csv")

# HXB
HXBers = c("HXB","BXH", "BN.","SHR")
get2 = inData_part2[,shortNames %in% HXBers]
selection = rowSums(get2 == "1/1") >= 1
out_HXBers = cbind(inData_part1, get2, inData_part3)[selection,]

#Write
write.csv(out_HXBers, "results/04_highImpact_HXB_BXH.csv")
out_HXBers = out_HXBers[!duplicated(paste(out_LEXFers$ID, out_LEXFers$TERM_NAME)),]
write.csv(out_HXBers, "results/04_highImpact_HXB_BXH_noDup.csv")


#LEXF
LEXFers = c("LEX","FXL","LXF","LE.","F34")
get2 = inData_part2[,shortNames %in% LEXFers]
selection = rowSums(get2 == "1/1") >= 1
out_LEXFers = cbind(inData_part1, get2, inData_part3)[selection,]

#Write
write.csv(out_LEXFers, "results/04_highImpact_LEXF_FXLE.csv")
out_LEXFers = out_LEXFers[!duplicated(paste(out_LEXFers$ID, out_LEXFers$TERM_NAME)),]
write.csv(out_LEXFers, "results/04_highImpact_LEXF_FXLE_noDup.csv")

#SS SR
Sers = c("SS.","SR.")
get2 = inData_part2[,shortNames %in% Sers]
selection = (rowSums(get2 == "1/1") >= 1) & (rowSums(get2 == "0/0") >= 1)
out_Sers = cbind(inData_part1, get2, inData_part3)[selection,]

#Write
write.csv(out_Sers, "results/04_highImpact_SS_SR.csv")
out_Sers = out_Sers[!duplicated(paste(out_Sers$ID, out_Sers$TERM_NAME)),]
write.csv(out_Sers, "results/04_highImpact_SS_SR_noDup.csv")

# LL LH LN
Lers = c("LL.","LH.","LN.")
get2 = inData_part2[,shortNames %in% Lers]
selection = (rowSums(get2 == "1/1") >= 1) & (rowSums(get2 == "0/0") >= 1)
out_Lers = cbind(inData_part1, get2, inData_part3)[selection,]

#Write
write.csv(out_Lers, "results/04_highImpact_LL_LN_LH.csv")
out_Lers = out_Lers[!duplicated(paste(out_Lers$ID, out_Lers$TERM_NAME)),]
write.csv(out_Lers, "results/04_highImpact_LL_LN_LH_noDup.csv")



## Cross compare research:
compare_me = read.csv("other_papers/LL_LN_LH_paper_results_compact.csv", header = T, sep = ";")
keep = compare_me[tolower(compare_me$Gene.name) %in% tolower(out_Lers$nearby_gene),]
keep2 = out_Lers[tolower(out_Lers$nearby_gene) %in% tolower(compare_me$Gene.name),] 


