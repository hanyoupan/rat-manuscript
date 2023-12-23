library(circlize)
library(data.table)
library(stringr)


#chromosomeList = paste("chr",c(1:20,"X"),sep = "")
chromosomeList = paste("chr",c(1:20),sep = "")

minQual = 0

############################################################################
#pdf("results/07_Circos_mRatBN7.pdf")

outnameSpace = 600000000

outNames = c("Deletions","Inversions","Duplications","Transpositions")
allTypes = c("DEL","INV","DUP", "TRA")
allColors = c("black","forestgreen","purple")
ratColors = c("darkgoldenrod1","darkgoldenrod3", "blue", "bisque3", "forestgreen","purple","red", "black")

circos_bar_height = 0.1
circos_bar_height_BN = 0.07
circos_text_sizes = 0.6
extraPlotThickness = 1



############################################################################
## Get ref sizes and positions

refGenomeSizes = read.csv("ref_genome_sizes.csv", sep = ";" ) 

allFactors = factor(c(chromosomeList, 
                    chromosomeList),
                    levels = chromosomeList)

allPositions = c(rep(0,length(chromosomeList)),
                 refGenomeSizes[1:length(chromosomeList),2])


############################################################################
## Initialize circos code
############################################################################

makeCircos = function(samples){

    #samples = samplesRN6
  
    first = TRUE
    circos.clear()
    circos.par("track.height" = 0.06, cell.padding = c(0.0, 0.4, 0.0, 0.4), start.degree = 90, canvas.xlim = c(-1.1,1.1), canvas.ylim = c(-1.2,1.2))
    circos.initialize(factors = allFactors, x = allPositions)
      
  
    for (i in 1:length(samples)){
        print(i)
        sample = samples[i]
        sampleColor = sampleColors[i]
        
        #################################################################################
        ## Read the VCF file
        #################################################################################
        print("Processing:")
        print(sample)
        
        
        circos.track(factors = allFactors, x = allPositions,
                     track.height = circos_bar_height_BN, ylim = c(-0.025,1.025),
                     panel.fun = function(x, y) {
                       if (i == 1){
                       circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(8, "mm"), 
                                   CELL_META$sector.index,facing = "clockwise", niceFacing = TRUE, 
                                   cex = 0.6, adj = c(1, 0.01))
                       }
                     })
        
        vcf1 <- fread(sample)
        
        
        ############################################################
        ## Quick test here to filter the results
        selectCol = colnames(vcf1)[ncol(vcf1)]
        vcf1 = vcf1[as.numeric(vcf1$QUAL) > minQual,]
        ############################################################
        

        colnames(vcf1)[1] = "chr" 
        
        chr1 = vcf1$chr
        vcf1 = vcf1[chr1 %in% chromosomeList]
        
        # check deletions only
        head(vcf1[vcf1$ALT != "<DEL>"])
        
        samples1 = colnames(vcf1)[10:ncol(vcf1)]
        

        
        #################################################################################
        # 01. Take deletions only
        vcf2 = vcf1#[vcf1$ALT == "<DEL>",]
        
        #################################################################################
        print(paste(nrow(vcf2), "total DELETIONS"))
        
        ####################################################
        # 02. SV type
        split4SV = strsplit(vcf2$INFO,"SVTYPE=")
        split4SV_2 = sapply(split4SV, "[", 2)
        split4SV_3 = sub(";.*","",split4SV_2)
        length(split4SV_3)
        
        SV_type = split4SV_3
        table(SV_type)
        vcf2$SV_type = SV_type
        
        ####################################################
        # 03. SV length
        split4SV = strsplit(vcf2$INFO,"SVLEN=")
        split4SV_2 = sapply(split4SV, "[", 2)
        split4SV_3 = sub(";.*","",split4SV_2)
        length(split4SV_3)
        
        SV_LEN = split4SV_3
        table(SV_type)
        vcf2$SV_LEN = SV_LEN
        
        ####################################################
        # 03. SV length / Version 1.
        split4SV = strsplit(vcf2$INFO,"END=")
        split4SV_2 = sapply(split4SV, "[", 2)
        split4SV_3 = sub(";.*","",split4SV_2)
        length(split4SV_3)
        
        SV_END = split4SV_3
        table(SV_type)
        vcf2$SV_END = as.numeric(SV_END)
        
        ####################################################
        # 03. SV length / Version 2.
        split4chrom = strsplit(vcf2$ALT,"\\[")
        split4chrom_2 = sapply(split4chrom, "[", 2)
        split4chrom_3 = sub(":.*","",split4chrom_2)
        chrom2 = split4chrom_3
        
        split4pos = strsplit(vcf2$ALT,":")
        split4pos_2 = sapply(split4pos, "[", 2)
        split4pos_3 = sub("\\].*|\\[.*","",split4pos_2)
        pos2 = split4pos_3
        pos2 = as.numeric(ifelse(is.na(pos2),0,pos2))
        
        # Q1. Are they all on the same chromosome?
        sum(chrom2 == vcf2$chr, na.rm = TRUE) # Yes.
        
        vcf2$SV_END = ifelse(is.na(SV_END), pos2, vcf2$SV_END)
        table(vcf2$SV_END  == pos2 | vcf2$ALT == "<DEL>")

        table(abs(as.numeric(vcf2$SV_LEN)) == abs((vcf2$POS - as.numeric(vcf2$SV_END))))
        
        vcf2$SV_LEN = abs((vcf2$POS - as.numeric(vcf2$SV_END)))
        
        sizeSV = sum(as.numeric(vcf2$SV_LEN))
        numberSV = nrow(vcf2)
        
        if (first == TRUE){
            outTable = rbind(c(numberSV, sizeSV))
            first = FALSE
        }else{
          outTable = rbind(outTable, rbind(c(numberSV, sizeSV)))
        }
        print(outTable)
        ####################################################
        # 06.Initialize Circos
        
        for(chr in unique(vcf2$chr)){
          draw_chr1 = vcf2$chr[vcf2$chr == chr]
          draw_pos1 = vcf2$POS[vcf2$chr == chr]
          draw_end1 = vcf2$SV_END[vcf2$chr == chr]
        
          #print(chr)
          #print(paste(length(draw_chr1), "DEL"))
          
          circos.rect(sector.index = chr,
                        xleft = draw_pos1,
                        ybottom = rep(0,length(draw_chr1)),
                        xright = draw_end1 + 2000000, #5000000 is pretty good
                        ytop = rep(1,length(draw_chr1)),
                        col = sampleColor,
                        border = NA)
        } 
        
    }    
    
    rownames(outTable) = sub("_rn6_large_svs.DELonly.vcf","",sub("data/","", samples))
    rownames(outTable) = sub("_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf","",rownames(outTable))
    colnames(outTable) = c("Number of large deletions","total length")
    return(outTable)
    
}  

###############################################################################
## Actually select samples and create circos.
###############################################################################

sampleColors = c("brown", "orange",rainbow(9))


samplesRN6 = paste0("data/",c(#"BN_male_rn6_large_svs.DELonly.vcf", # 0
                            "BN-Lx_Cub_rn6_large_svs.DELonly.vcf",
                            "SHR_OlaIpcv_rn6_large_svs.DELonly.vcf", #1
                            "BXH10_rn6_large_svs.DELonly.vcf",
                            "BXH8_rn6_large_svs.DELonly.vcf",# 3
                            "BXH2_rn6_large_svs.DELonly.vcf", 
                            "HXB17_rn6_large_svs.DELonly.vcf",# 5
                            "HXB2_rn6_large_svs.DELonly.vcf", 
                            "HXB21_rn6_large_svs.DELonly.vcf" # 7
                            
))
samplesBN7 = paste0("data/",c(#"BN_male_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf",# 0
                              "BN-Lx_Cub_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf",
                            "SHR_OlaIpcv_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf", #1
                            "BXH10_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf", 
                            "BXH8_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf", # 3
                            "BXH2_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf", 
                            "HXB17_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf", # 5
                            "HXB2_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf",
                            "HXB21_mRatNor1_mRatBN7.2_large_svs.DELonly.vcf" # 7
                            
))

pdf(paste("results/01_circos_plots_Q",minQual,".pdf", sep = ""))
# Total number of large DELs
RN6  = makeCircos(samplesRN6) # RN6
BN72 = makeCircos(samplesBN7) # BN7.2
dev.off()

write.csv(RN6, paste("results/RN6_numberOfDeletions_Q",minQual,".csv",sep = ""))
write.csv(BN72, paste("results/BN72_numberOfDeletions_Q",minQual,".csv",sep = ""))

# show the tables
RN6
BN72


