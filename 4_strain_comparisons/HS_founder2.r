
library(UpSetR)

chrList = paste("chr",c(1:20,"X","Y"),sep = "")
#chrList = paste("chr",c(1:5),sep = "")

teal = "#31A9B8"
salmon = "#F8766D"

# source: https://www.ncbi.nlm.nih.gov/grc/rat/data
chrSizes = c(260522016  ,
              249053267 ,
              169034231 ,
              182687754 ,
              166875058 ,
              140994061	,
              135012528 ,
              123900184	,
              114175309	,
              107211142	,
              86241447	,
              46669029	,
              106807694	,
              104886043	,
              101769107	,
              84729064	,
              86533673	,
              83828827	,
              57337602	,
              54435887	,
              152453651	,
              18315841)



dataOrder = c("ACI.N_JL.ILM","ACI.N_AP.ILM",
              "BN.N_JL.ILM","BN.N_AP.ILM",
              "BUF.N_JL.ILM","BUF.N_AP.ILM",
              "F344.N_JL.ILM","F344.N_AP.ILM",
              "M520.N_JL.ILM","M520.N_AP.ILM",
              "MR.N_JL.ILM","MR.N_AP.ILM",
              "WN.N_JL.ILM","WKY.N_AP.ILM",
              "WKY.N_JL.ILM","WN.N_AP.ILM")

strains = substr(dataOrder,1,3)

first = TRUE
for (chr in chrList){
    print(paste("processing:", chr))
    inData = read.csv(paste("results/foundInAll_",chr,".gvcf.csv.founder_analysis",sep = ""), sep = "\t", header = T)
    
    inHigh = inData[inData$impact == "HIGH",]
  
    # All calls 
    inCall_information = inData[,1:8]
    inCall = inData[,9:(ncol(inData)-5)]
    inCall = inCall[,dataOrder]
    inOther = inData[,(ncol(inData)-4):ncol(inData)]

    # High impact data information
    inCall_High = inHigh[,9:(ncol(inData)-5)]  
    inOther_High = inHigh[,(ncol(inData)-4):ncol(inData)]

    # Convert to binary
    inCall_binary = inCall
    filter = rowSums(inCall_binary == "1/1") >= 1
    
    table(filter)
    inCall_binary = inCall_binary[filter,]
    inCall_information = inCall_information[filter,]
    inOther =  inOther[filter,]

    # Fold the data into a single table... if it fits.
    if (first == TRUE){
      outData_gt = inCall_binary
      outData_information = inCall_information
      outData_other = inOther
      first = FALSE
    }else{
      outData_gt = rbind(outData_gt,  inCall_binary)
      outData_information = rbind(outData_information, inCall_information)
      outData_other = rbind(outData_other, inOther)
    }
}


outData_binary = outData_gt == "1/1"
outData_unique = outData_binary[rowSums(outData_binary) <=2,]


##################################################
## Prepare data for new figure
##################################################

# convert to binary
outData_gt_binary = outData_gt == "1/1"

# Merge per strain
strains = substr(colnames(outData_unique),1,3)

unique(strains)

mergeStrains = function(strain){
  outNumbers = rowSums(outData_gt_binary[,strains==strain])
  return(outNumbers)
}

ACI = mergeStrains("ACI")
BN = mergeStrains("BN.")
BUF = mergeStrains("BUF")
F344 = mergeStrains("F34")
M520 = mergeStrains("M52")
MR = mergeStrains("MR.")
WN = mergeStrains("WN.")
WKY = mergeStrains("WKY")

outData_strains = cbind(ACI, BN,BUF,F344,M520,MR,WN,WKY)

ACI = 0
BN = 0
BUF = 0
F344 = 0
M520 = 0
MR = 0
WN = 0
WKY = 0

## Unique to strain
filter = rowSums(outData_strains >= 1) == 1
outData_strains_u = outData_strains[filter,]
outData_information_u = outData_information[filter,]

filter = rowSums(outData_strains >= 1) == 8
outData_strains_a = outData_strains[filter,]
outData_information_a = outData_information[filter,]

nrow(outData_strains)
nrow(outData_information)

##################################################
## Create plot 1. Number of SNPs per strain
##################################################
pdf("founders_figure1.pdf", width = 10, height = 10)
par(mar = c(5.1, 8.1, 4.1, 2.1) ) # c(bottom, left, top, right))

# BN, MR, BUF, WN, M520, F344, ACI, WKY
plotColors = c("brown","lightblue","#651b20","orange","pink","ivory","#5d594f","#c8ac97")

figure1 = colSums(outData_strains >= 1)
figure1b = colSums(outData_strains_u >= 1)

shared = table(figure1c)

plotOrder = order(figure1)
barplot((figure1[plotOrder])/1000000,
        ylab = "Total number of SNPs (millions)\n", 
        xlab = "",
        main = "Total number of SNPs per HRDP founder",
        col = plotColors,
        las = 1,
        cex.axis = 1.2,
        cex.names = 1.2)

bar = barplot((figure1b[plotOrder])/1000000, add = T, col =   adjustcolor( "grey", alpha.f = 0.5), las = 1)
text(x = bar[2:length(bar)], y = rep(0.1,7), "unique", cex = 0.7)

#######################################################
## Alternate version
#######################################################

plotOrder = order(figure1)
barplot((figure1[plotOrder])/1000000,
        yaxt = "n",
        ylab = "Total number of SNPs\n", 
        xlab = "",
        main = "Total number of SNPs per HRDP founder",
        col = plotColors,
        cex.axis = 1.4,
        cex.names = 1.4,
        las = 1)

axis(side = 2, at = seq(0,5,1), cex.axis = 1.4, labels = paste(seq(0,5,1), "M"), las = 2)

bar = barplot((figure1b[plotOrder])/1000000, add = T, xaxt = "n", yaxt = "n", col =   adjustcolor( "grey", alpha.f = 0.5), las = 1)
text(x = bar[2:length(bar)], y = rep(0.1,7), "unique", cex = 1)



dev.off()
##################################################
# plot 2.
##################################################
pdf("founders_figure2.pdf")
par(mar = c(5.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right))

plotColors = c("#9fe3fe", "#9cb5fe",
               "#df9bfd", "#fd9dcb", 
               "#fd9dcb", "#fecf9f",
               "#fee9a7", "#f7fca8",
               "#d3ffa5")

unique(strains)

outData_shared = rowSums(outData_strains >= 1)

# try shared with at least
plotTable0 = table(outData_information$X.CHROM)
b = barplot(plotTable0[chrList]/1000, las = 2, add = F,col = "grey", yaxt = 'n', 
            cex.axis = 1.2, cex.names = 1.2,
            ylab = "Total number of SNPs\n")
axis(side = 2, at = seq(0,1400,200), labels = paste(seq(0,1.400,.200), "M"), cex.axis = 1.2,  las = 2)
axis(side = 2, at = seq(0,1400,200), labels = paste(seq(0,1400,200), "M"), cex.axis = 1.2,  las = 2)



for (i in 1:length(unique(strains))){
  plotTable1 = table(outData_information$X.CHROM[outData_shared >= i])
  barplot(plotTable1[chrList] /1000, las = 2, add = T,col = plotColors[i], yaxt = "n", xaxt = "n")
   print(i)
}

legend(x = "topright",
       legend = paste(c(1,2,3,4,5,6,7,8), c("strain",rep("strains",7)),sep = " "),
       title = "Found in:",
       col = plotColors, pch = 15, bty = "n", cex = 1.4)

#########################################################
## Different labelling method
#########################################################

# try shared with at least
plotTable0 = table(outData_information$X.CHROM)
b = barplot(plotTable0[chrList]/1000, las = 2, add = F,col = "grey", ylab = "total number of SNPs (thousands)")
#axis(side = 2, at = seq(100,1200,100), labels = paste(seq(100,1200,100),), las = 2)


for (i in 1:length(unique(strains))){
  plotTable1 = table(outData_information$X.CHROM[outData_shared >= i])
  barplot(plotTable1[chrList] /1000, las = 2, add = T,col = plotColors[i], xaxt = "n", yaxt = "n")
  print(i)
}

legend(x = "topright",
       legend = paste(c(1,2,3,4,5,6,7,8), c("strain",rep("strains",7)),sep = " "),
       title = "Found in:",
       col = plotColors, pch = 15, bty = "n", cex = 1.2)

outData_strains_u
outData_information_u

dev.off()

##################################################
# plot 3.
##################################################

filter = rowSums(outData_strains >= 1) == 1 | rowSums(outData_strains >= 8)

outData_final = outData_strains[filter,] >= 1 

uniqueToStrains = colSums(outData_strains_u >= 1)
sharedByStrains = colSums(outData_strains_a >= 1)

head(outData_strains_a)
head(outData_information_a)


write.csv(cbind(outData_information_a, outData_strains_a), "sharedByFounders.csv")


##################################################
# Number of high impact SNPs
##################################################

head(outData_gt)

head(outData_information)

table(outData_information$impact)

head( outData_other)

max(outData_other$all_alt)


