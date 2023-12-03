

## 00. Load Packages
library(ggplot2)


chrList = paste("chr",c(1:20, "X","Y"),sep = "")
#chrList = paste("chr",c(1:2),sep = "")

pdf(file = "results/01_HXB_checkup.pdf", width =10, height = 5)

shrSamp = c(40,41)
bnSamp = c(1,2)
setOrder = c( 1,2, 29,  8,23, 22,  7, 14,  3,  5, 39, 13, 12,  6, 30, 20, 18, 11, 26, 21, 16, 25, 28, 24,19,38,35,10,34,15,32,9,27,17,31,33,4,37,36,40,41)

grandTable = c()
First = TRUE

doKaryoPlots = TRUE

barplotColors = c("red", "yellow", "forestgreen","blue","grey", "magenta")



################################################################################
## 01. Load in 1 chromosome.

for(chr in chrList){
    print(paste("Processing", chr))
    #chr = "chr1"
    currentChromosome = chr
    
    quality = read.csv(paste("data/deepvariant_",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HXBonly.gvcf.gq",sep = ""), header = T, sep = "\t", row.names =1)
    #genotype = read.csv(paste("data/deepvariant_",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HXBonly.gvcf.gt",sep = ""), header = T, sep = "\t", row.names =1)
    
    #quality = read.csv(paste("data/deepvariant_",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HXBonly.gvcf.gq",sep = ""), header = T, sep = "\t", row.names =1)
    genotype = read.csv(paste("data/deepvariant_",chr,"_163_samples_mRatBN7.2_annot_error_removed_GQ30_HXBonly.gvcf.adj",sep = ""), header = T, sep = "\t", row.names =1)
    
    coordinate = sub("chr.._*", "", rownames(genotype))   
    coordinate = sub("_.*", "", coordinate)  
    genotype$coordinates = as.numeric(coordinate)
    
    head(quality)
    head(genotype)
    
    shortNames = colnames(quality)
    ################################################################################
    ## 02. Identify from which parent the call came.
    head(genotype)
    
    # Get which calls are homozygous in either BN or SHR
    select_BN_1 = rowSums(genotype[,bnSamp] == "1/1") >= 1
    select_SHR_1 = rowSums(genotype[,shrSamp] == "1/1") >= 1
    
    # Also select where it has to be in both
    select_BN_2 = rowSums(genotype[,bnSamp] == "1/1") >= 2
    select_SHR_2 = c(rowSums(genotype[,shrSamp] == "1/1") >= 2)
    
    # See where the HXB calls get their SNPs
    HXB_BNsource = genotype[select_BN_1 & !select_SHR_1,]
    HXB_SHRsource = genotype[select_SHR_1 & !select_BN_1,]
    HXB_BOTHsource = genotype[select_SHR_1 & select_BN_1,]
    HXB_DENOVOsource = genotype[!select_SHR_1 & !select_BN_1,] 
    
    ## Select the homozygous SNPs
    colSums(HXB_BNsource == "1/1")
    colSums(HXB_SHRsource == "1/1")
    colSums(HXB_SHRsource == "0/1")
    colSums(HXB_BOTHsource == "1/1")
    colSums(HXB_DENOVOsource == "1/1")
    
    
    ################################################################################
    ## 03. Find a way to create an overview, then a summation
    
    overview_1 = cbind((colSums(HXB_BNsource == "1/1")),
                       (colSums(HXB_SHRsource == "0/0")),
                       (colSums(HXB_SHRsource == "0/1")),
                        (colSums(HXB_SHRsource == "1/1")),
                        (colSums(HXB_BOTHsource == "1/1")),
                        ( colSums(HXB_DENOVOsource == "1/1")))
    colnames(overview_1) = c("BN novo","BN", "HET","SHR","Both","DeNovo")
    
    # Total number of SNPs by source, relative to the reference

    
    if (First == TRUE){
      grandTable = overview_1
      #setOrder = order(overview_1[,4])
    }else{
      grandTable = grandTable + overview_1
    }
    overview_1 = overview_1[setOrder,]
    
    ################################################################################
    ## 03. Plot everything for visual thinkers:
    

    par(mar = c(7.1, 4.1, 4.1, 2.1)) # c(bottom, left, top, right)
    
    # Barplot for visualisation
    barplot(as.matrix(t(overview_1)/1000), las = 2, 
            col = barplotColors, 
            names = sub("_.*", "",rownames(overview_1)), 
            cex.names = 1,
            main = currentChromosome,
            ylab = "total number of SNPs x1000", ylim = c(0,max(as.matrix(t(overview_1)/1000)*1.5)))
    
    legend(x = "top", title = "Source of SNP:", legend = colnames(overview_1), pch = 15, 
           col = barplotColors, horiz = T, cex = 1, bty = "n")
    
    
    ################################################################################
    ## 03b. Plot everything per position
    
    if (doKaryoPlots == TRUE){
        print("making position dependent plot...")
        #coordinate = sub("chr.._*", "", rownames(HXB_SHRsource))   
        #coordinate = sub("_.*", "", coordinate)   
        #coordinate = as.integer(coordinate)
        
        HXB_SHRsource$coordinates

        pdf(file = paste0("results/karyogram",chr,".pdf"), width =20, height = 12)
        
        howManyPlot = 41
        barBroad = 1.2
        shortNames = sub("_.*", "",colnames(HXB_SHRsource))
        
        ## Create plotting region
        plot(x = NA,xlim = c(0,(howManyPlot+1)*barBroad), ylim = c(0,max(genotype$coordinates)* 1.1), xaxt='n', main = chr)
        Axis(side=1, at = 1:howManyPlot*barBroad, labels=shortNames[1:howManyPlot], las = 2, cex.axis=0.5)
        
        legend(x = "top", title = "Source of SNP:", legend = colnames(overview_1), pch = 15, 
               col = barplotColors, horiz = T, cex = 1, bty = "n")
        
        
        barWidth = 0.15
        ## Function for easier drawing
        drawForMe = function(inCoords, inColor, inOffset, fatness ){
          
          #inCoords = drawCoordinates1 
          #inColor = barplotColors[1] 
          #inOffset = -0.5
          
          rect(xleft = rep(i*barBroad,length(inCoords)) + inOffset + 0, 
               ybottom = inCoords, 
               xright = rep(i*barBroad,length(inCoords)) +  inOffset + barWidth, 
               ytop = inCoords + fatness, 
               col = inColor, border = NA)
        }
        
        ## Call the function for every sample and type of SNP, like an amateur
        for (i in 1:howManyPlot){
          
          drawCoordinates1 = HXB_BNsource$coordinates[HXB_BNsource[,i] == "1/1"]
          drawCoordinates2 = HXB_SHRsource$coordinates[HXB_SHRsource[,i] == "0/0"]
          drawCoordinates3 = HXB_SHRsource$coordinates[HXB_SHRsource[,i] == "0/1"]
          drawCoordinates4 = HXB_SHRsource$coordinates[HXB_SHRsource[,i] == "1/1"]
          #drawCoordinates5 = HXB_BOTHsource$coordinates[HXB_BOTHsource[,i] == "1/1"]
          drawCoordinates6 = HXB_DENOVOsource$coordinates[HXB_DENOVOsource[,i] == "1/1"]
          

          drawForMe(inCoords = drawCoordinates1, inColor = barplotColors[1], inOffset = -0.45, fatness =  100)
          drawForMe(inCoords = drawCoordinates2, inColor = barplotColors[2], inOffset = -0.30, fatness =  100)
          drawForMe(inCoords = drawCoordinates3, inColor = barplotColors[3], inOffset = -0.15, fatness =  100)
          drawForMe(inCoords = drawCoordinates4, inColor = barplotColors[4], inOffset = 0.0, fatness =  100)
          #drawForMe(inCoords = drawCoordinates5, inColor = barplotColors[5], inOffset = 0.15)
          drawForMe(inCoords = drawCoordinates6, inColor = barplotColors[6], inOffset = 0.15, fatness =  100)

        }
        dev.off()
    
    }
  ################################################################################
    
    First = FALSE
}

################################################################################
## 04. Establish the true order
trueOrder = order(grandTable[,4])
#grandTable = grandTable[order(grandTable[,4])]




################################################################################
## 05. Write the total out table.
grandTable = grandTable[trueOrder,]
grandTable = grandTable[-3,] # Remove the coordinates line

write.csv(grandTable,"results/02_SNP_origin_HXB_SHRorBN.csv")

barbar = barplot(as.matrix(t(grandTable)/1000), las = 2, 
        col = barplotColors, 
        names = sub("_.*", "",rownames(grandTable)), 
        cex.names = 1,
        main = "All chromosomes",
        ylab = "total number of SNPs x1000", ylim = c(0,max(as.matrix(t(grandTable)/1000)*1.4)))

legend(x = "top", title = "Source of SNP:", legend = colnames(grandTable), pch = 15, 
       col = barplotColors, horiz = T, cex = 1, bty = "n")

trueOrder = row.names(grandTable)

dev.off()


## Do a little update to the table
grandTable2 = grandTable
colnames(grandTable2) = paste(colnames(grandTable2), c("(1/1)","(0/0)","(0/1)","(1/1)","(1/1)","(1/1)"))
write.csv(grandTable2,"results/02_SNP_origin_HXB_SHRorBN_v2.csv")





