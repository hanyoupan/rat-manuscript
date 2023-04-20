library("rio")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")
library("gridExtra")
library("ggpubr")


# change csv files to df, preprocess them and clean rfid columns 
csv_to_df <- function(csv_file){
  df_temp<-read.table(file=csv_file, header=T, quote="", sep=",")
  return(df_temp)
}

# read the files in rawdata folder, convert them to csv
rawfile_list = list.files("./all_summary", pattern=".csv", all.files=FALSE, full.names=FALSE)

rawfile_list[1]

add_ref_name<-function(sample_name){
  if(grepl('_summary', sample_name,fixed=T)){
    split_file_names <-strsplit(sample_name,'_summary')[[1]][1]
    ref <- split_file_names %>% strsplit('_') %>% sapply(tail, 1)
    ref
  }
  else{
    split_file_names <-strsplit(sample_name,'summary')[[1]][1]
    ref <- split_file_names %>% strsplit('_') %>% sapply(tail, 1)
    ref
  }
}

add_strain_name<-function(sample_name){
  if(grepl('rn6', sample_name,fixed=T)){
    #split_file_names <-strsplit(sample_name,'rn6')[[1]][1]
    strain <- sample_name %>% strsplit('_rn6') %>% sapply(head, 1)
    strain
  }
  else{
    #split_file_names <-strsplit(sample_name,'mRatNor1')[[1]][1]
    strain <- sample_name %>% strsplit('_mRatNor1') %>% sapply(head, 1)
    strain
  }
}
   
df0 <- data.frame(matrix(ncol = 33, nrow = 0))
for(i in 1:length(rawfile_list)){ 
  cat(rawfile_list[i])
  #print('/n')
  #cat(add_ref_name(rawfile_list[i]))
  print('/n')
  file_name = paste0("./all_summary", "/", rawfile_list[i])
  df_temp <- csv_to_df(file_name)
  df_temp$ref <- add_ref_name(rawfile_list[i])
  df_temp$strain <- add_strain_name(rawfile_list[i])
  df0 <- rbind(df_temp[,1:33],df0[,1:33])
}

df0 <- df0 %>% select(ref,strain,everything())

df0$ref <- gsub('mRatNor1','mRatBN7.2',df0$ref)

df0$ref <- gsub('rn6','Rnor_6.0',df0$ref)
 

dim(df0)  

write.csv(df0, "all.csv")

#df0_strain <- mutate(df0, sitevar=paste(ref,strain,sep=", "))

dif_rn6_mRatBN7.2 <- setdiff(subset(df0,(df0$ref=='Rnor_6.0'))$strain,unique(subset(df0, df0$ref=='mRatBN7.2')$strain))
dif_mRatBN7.2_rn6 <- setdiff(unique(subset(df0, df0$ref=='mRatBN7.2')$strain),subset(df0,(df0$ref=='Rnor_6.0'))$strain)

df0_both <- subset(df0, !((df0$strain %in% dif_rn6_mRatBN7.2)| (df0$strain %in% dif_mRatBN7.2_rn6) ))

#write.csv(df0_both, "df_both.csv")

col_var <- colnames(df0_both)

pdf(file="comparison_point_graphs.pdf", width=11, height=8.5, paper="special")
for(i in (5:length(colnames(df0_both)))){
  plot<-ggplot(df0_both) + geom_point(aes_string(x=col_var[i], y=paste("reorder(",col_var[2],",",col_var[i],")"), color=col_var[1], size = 3))+ylab("strain")
  plot = plot + guides(size = FALSE) + theme_bw(base_size=14)+ guides(colour = guide_legend(override.aes = list(size = 5)))+theme(legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.key = element_rect(fill=NA))
  
  print(plot)
  #ggplot(df0_both) + geom_point(aes(x=gems_detected, y=reorder(strain,gems_detected), color=ref))
#  ggplot(df0_both, aes(x=reorder(ref,gems_detected), y=gems_detected, fill=ref)) + geom_boxplot(fill="white")+ geom_dotplot(binaxis='y', stackdir='center',binwidth = 25000) + xlab("ref")
}
dev.off()
       


rn6<-read.table(file="./deepvariant_rn6_36_samples.cnt", header=T, sep="\t")
bn7<-read.table(file="./deepvariant_mRatNor1_36_samples.cnt", header=T, sep="\t")
#save(file="deepvar_rn6_bn7_36samples.Rdata", rn6, bn7) 

load("deepvar_rn6_bn7_36samples.Rdata")



rn6$ref<-"Rnor_6.0"
bn7$ref<-"mRatBN7.2"
df0<-rbind(rn6,bn7)

# mis==Missing
df0<-subset(df0, Qual>=30)

#chrs<-c(paste0("chr",1:20) , "chrX", "chrY")
#df0$chr<- factor(df0$chr,levels = chrs) 
#df0<-subset(df0, chr!="chrX" )
#df0<-subset(df0, chr!="chrY")
#chrs<-c( paste0("chr",1:20))
#df0$chr<- factor(df0$chr,levels = chrs) 
#table(df0$chr)
snp<-subset(df0,Type=="SNP")
indel<-subset(df0,Type=="Indel")
dim(snp)
dim(indel)

s0=10
s1 =12
head(snp)
#altsnp <- ggplot(subset(snp), aes(x=chr, y=alt, fill=ref)) + geom_violin(scale="width") + labs(y="Count", x="", title = "Homozygotic alternative SNPs from 36 inbred strains")+ theme(axis.text.x = element_text(size=s0, angle=45))
#altindel <- ggplot(subset(indel), aes(x=chr, y=alt, fill=ref)) + geom_violin(scale="width") + labs(y="Count", x="", title = "Homozygotic alternative Indels from 36 inbred strains") +theme(axis.text.x = element_text( size =s0, angle=45))

#entire genome
names(snp)
altsnpvio0   <- ggplot(snp, aes(x=ref, y=AltCnt, fill=ref))   + theme_bw() + geom_violin(adjust=3, scale="area", show.legend=F) + labs(y="Count", x="", title = "Homozygous SNPs")+ ylab("Strains")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1))
altindelvio0 <- ggplot(indel, aes(x=ref, y=AltCnt, fill=ref)) + theme_bw() + geom_violin(adjust=2, scale="area", show.legend=F) + labs(y="Count", x="", title = "Homozygous Indels") +ylab("Strains") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1))
HETsnpvio0   <- ggplot(snp, aes(x=ref, y=HetCnt, fill=ref))   + theme_bw() + geom_violin(adjust=3, scale="area", show.legend=F) + labs(y="Count", x="", title = "Heterozygous SNPs")+ ylab("Strains")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1))
HETindelvio0 <- ggplot(indel, aes(x=ref, y=HetCnt, fill=ref)) + theme_bw() + geom_violin(adjust=2, scale="area", show.legend=F) + labs(y="Count", x="", title = "Heterozygous Indels") +ylab("Strains") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1))

## FLIPPED
altsnpvio0   <- ggplot(snp, aes(x=ref, y=AltCnt, fill=ref))   + theme_bw() + geom_violin(adjust=3, scale="area", show.legend=F) + labs(y="Count", x="", title = "Homozygous SNPs")+ ylab("Strains")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))
altindelvio0 <- ggplot(indel, aes(x=ref, y=AltCnt, fill=ref)) + theme_bw() + geom_violin(adjust=2, scale="area", show.legend=F) + labs(y="Count", x="", title = "Homozygous Indels") +ylab("Strains") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))
HETsnpvio0   <- ggplot(snp, aes(x=ref, y=HetCnt, fill=ref))   + theme_bw() + geom_violin(adjust=3, scale="area", show.legend=F) + labs(y="Count", x="", title = "Heterozygous SNPs")+ ylab("Strains")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))
HETindelvio0 <- ggplot(indel, aes(x=ref, y=HetCnt, fill=ref)) + theme_bw() + geom_violin(adjust=2, scale="area", show.legend=F) + labs(y="Count", x="", title = "Heterozygous Indels") +ylab("Strains") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))


head(snp)

#where is alt coming from?

min_strain=36
# comAltSNP<-aggregate(AltCnt~ref+chr, FUN=length, data=subset(snp, AltCnt>=min_strain))
# comAltIndel<-aggregate(AltCnt~ref+chr, FUN=length, data=subset(indel, AltCnt>=min_strain ))
# comHETSNP<-aggregate(HetCnt~ref+chr, FUN=length, data=subset(snp, HetCnt>=min_strain ))
# comHETIndel<-aggregate(HetCnt~ref+chr, FUN=length, data=subset(indel, HetCnt>=min_strain ))

length1M<-function(x){length(x)/1e6}
length1K<-function(x){length(x)/1000}

names(snp)
# entire genome, SNPs/Indels found in all 36 strains
comAltSNP0<-aggregate(AltCnt~ref, FUN=length1K, data=subset(snp, AltCnt>=min_strain ))
comAltIndel0<-aggregate(AltCnt~ref, FUN=length1K, data=subset(indel, AltCnt>=min_strain ))
comHETSNP0<-aggregate(HetCnt~ref, FUN=length1K, data=subset(snp, HetCnt>=min_strain ))
comHETIndel0<-aggregate(HetCnt~ref, FUN=length1K, data=subset(indel, HetCnt>=min_strain ))
#comMisSNP0<-aggregate(MissCnt~ref, FUN=length, data=subset(snp, MissCnt>=min_strain ))
#comMisIndel0<-aggregate(MissCnt~ref, FUN=length, data=subset(indel, MissCnt>=min_strain))

# number, in millions, of variants per reference genome
AltSNP0<-aggregate(AltCnt~ref, FUN=length1M, data=subset(snp, AltCnt>0))
AltIndel0<-aggregate(AltCnt~ref, FUN=length1M, data=subset(indel, AltCnt>0))
HETSNP0<-aggregate(HetCnt~ref, FUN=length1M, data=subset(snp, HetCnt>0))
HETIndel0<-aggregate(HetCnt~ref, FUN=length1M, data=subset(indel, HetCnt>0))

#HOM SNPs, 36 strains
p.AltSNP0<-ggplot(AltSNP0, aes(x=ref, y=AltCnt, fill=ref))      + theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous SNP")    +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1)) + coord_cartesian(ylim=c(0,10))   + ylab("Count, x 1M")
p.HETSNP0<-ggplot(HETSNP0, aes(x=ref, y=HetCnt, fill=ref))      + theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous SNP")    +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1)) + coord_cartesian(ylim = c(0,5)) + ylab("Count, x 1M")
p.AltIndel0<-ggplot(AltIndel0, aes(x=ref, y=AltCnt, fill=ref))  + theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous Indels") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1)) + coord_cartesian(ylim = c(0,10)) + ylab("Count, x 1M")
p.HETIndel0<-ggplot(HETIndel0, aes(x=ref, y=HetCnt, fill=ref))  + theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous Indels") +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0), plot.title = element_text(size=s1)) + coord_cartesian(ylim = c(0,5)) + ylab("Count, x 1M")

p.comAltSNP0<-ggplot(comAltSNP0, aes(x=ref, y=AltCnt, fill=ref))       + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous SNPs, all 36 strains")        +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comAltIndel0<-ggplot(comAltIndel0, aes(x=ref, y=AltCnt, fill=ref))   + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous Indels, all 36 strains")  +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comHETSNP0<-ggplot(comHETSNP0, aes(x=ref, y=HetCnt, fill=ref))       + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous SNPs,  all 36  strains")      +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comHETIndel0<-ggplot(comHETIndel0, aes(x=ref, y=HetCnt, fill=ref))   + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous Indels, all 36 strains") +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0), plot.title = element_text(size=s0)) + ylab("Count, x 1K")


#flipped coords
p.AltSNP0<-ggplot(AltSNP0, aes(x=ref, y=AltCnt, fill=ref))      + ylim(c(0,10))+ theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous SNP")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))+ ylab("Count, x 1M")
p.HETSNP0<-ggplot(HETSNP0, aes(x=ref, y=HetCnt, fill=ref))      + ylim(c(0,5))+ theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous SNP")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1)) +ylab("Count, x 1M")
p.AltIndel0<-ggplot(AltIndel0, aes(x=ref, y=AltCnt, fill=ref))  + ylim(c(0,10))+ theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous Indels")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1)) + ylab ("Count, x 1M")
p.HETIndel0<-ggplot(HETIndel0, aes(x=ref, y=HetCnt, fill=ref))  + ylim(c(0,5))+ theme_bw() + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous Indels")  +theme(axis.text.x = element_text(size =s0),axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s1))+ylab("Count, x 1M")

p.comAltSNP0<-ggplot(comAltSNP0, aes(x=ref, y=AltCnt, fill=ref))      + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous SNPs, all 36 strains")        +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comAltIndel0<-ggplot(comAltIndel0, aes(x=ref, y=AltCnt, fill=ref))  + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Homozygous Indels, all 36 strains")  +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comHETSNP0<-ggplot(comHETSNP0, aes(x=ref, y=HetCnt, fill=ref))      + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous SNPs,  all 36  strains")      +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s0))+ylab("Count, x 1K")
p.comHETIndel0<-ggplot(comHETIndel0, aes(x=ref, y=HetCnt, fill=ref))  + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Heterozygous Indels, all 36 strains") +theme(axis.text.x = element_text(size =s0), axis.text.y = element_text(size =s0, angle = 90, hjust=0.5), plot.title = element_text(size=s0)) + ylab("Count, x 1K")

#p.comMisSNP0<-ggplot(comMisSNP0, aes(x=ref, y=mis, fill=ref)) + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="SNP not called in at least 18 strains", size=1)  +theme(axis.text.x = element_text(size =s0), plot.title = element_text(size=12))
#p.comMisIndel0<-ggplot(comMisIndel0, aes(x=ref, y=mis, fill=ref)) + geom_bar(stat="identity", position="dodge", show.legend=F) + labs(y="Count", x="", title="Indel not called in at least 18 strains")  +theme(axis.text.x = element_text(size =s0), plot.title = element_text(size=12))

## Generate figure 1. part 1.##  +theme(axis.text.x = element_text(size =s0), plot.title = element_text(size=12))
pointSize = 0.2


stored_names = df0_both[,"strain"]
df0_both[,"strain"][c(1,2)] = c("SHR1","SHR1")
df0_both[,"strain"][c(63,64)] = c("BN1","BN1") ## BN_male x2
df0_both[,"strain"][c(65,66,67,68,69,70,71,72)] = c("BN2","BN2","BN3","BN3","BN4","BN4","BN5","BN5") # BN-LX, 0503_BN-male, 0502_BN-eva, 0502_BN-eva-clif

AdjustedNames = cbind(stored_names,df0_both[,"strain"])

i = 15+4
plot<-ggplot(df0_both, size=0.5) + geom_point(aes_string(x=col_var[i], y=paste("reorder(",col_var[2],",",col_var[i],")"), color=col_var[1])) + labs(y="", x="Mean Depth", title = "Mean Depth of Coverage,  36 strains")+ ylab("")
plot1 = plot + guides(size = FALSE) + coord_flip() + theme_bw()+ guides(colour = guide_legend(override.aes = list(size = 5)))+ theme(legend.position = "none", legend.title = element_text(size = s0),  plot.title = element_text(size=s1), axis.title.y=element_text(size=s0), axis.text.x = element_text(angle = 90, vjust=0.5), legend.text = element_text(size = s0), legend.key = element_rect(fill=NA))
i = 16+4
plot<-ggplot(df0_both, size=0.5) + geom_point(aes_string(x=col_var[i], y=paste("reorder(",col_var[2],",",col_var[i],")"), color=col_var[1])) + labs(y="", x="% Zero coverage", title = "Zero Coverage,  36 strains")+ ylab("")
plot2 = plot + guides(size = FALSE) + coord_flip() + theme_bw()+ guides(colour = guide_legend(override.aes = list(size = 5)))+ theme(legend.position = "none", legend.title = element_text(size = s0), plot.title = element_text(size=s1),  axis.title.y=element_text(size=s0), axis.text.x = element_text(angle = 90, vjust=0.5), legend.text = element_text(size = s0), legend.key = element_rect(fill=NA))
i = 17+4
plot<-ggplot(df0_both, size=0.5) + geom_point(aes_string(x=col_var[i], y=paste("reorder(",col_var[2],",",col_var[i],")"), color=col_var[1]))+ labs(y="", x="% Mapped Reads", title = "Mapped Reads,  36 strains")+ ylab("")
plot3 = plot + guides(size = FALSE) + coord_flip() + theme_bw()+ guides(colour = guide_legend(override.aes = list(size = 5)))+ theme(legend.position = "none", legend.title = element_text(size = s0), plot.title = element_text(size=s1),  axis.title.y=element_text(size=s0),axis.text.x = element_text(angle = 90, vjust=0.5), legend.text = element_text(size = s0), legend.key = element_rect(fill=NA))
  

## This is where the figure for the paper is created.

## My blue
myBlue = "#00bfc4"
myRed = "#f8766d"

#df0_both$ref <- factor(df0_both$ref,
#                       levels = c("Rnor_6.0",'mRatBN7.2'),ordered = TRUE)

## Step 1. Turn plot2 and plot3 into boxplots... somehow
new2 <- ggplot(df0_both, aes(x=ref, y=zero_coverage)) + theme_bw() + 
        labs(y="Zero coverage, %", x="", title="Zero coverage, %")  + 
        theme(legend.position = "none", legend.title = element_text(size = s1), 
         plot.title = element_text(size=s1),  axis.title.y=element_text(size=s1),
         axis.text.x = element_text(angle = 0, vjust=0.5),
         axis.text.y = element_text(angle = 90, vjust=0.05, hjust = 0.5),legend.text = element_text(size = s1), 
         legend.key = element_rect(fill=NA)) +
         geom_boxplot(color="black", fill=c(myBlue,myRed))

  
new2

new3 <- ggplot(df0_both, aes(x=ref, y=mapped_reads)) + theme_bw() + 
  labs(y="Mapped reads, %", x="", title="Mapped reads, %") +
  theme(legend.position = "none", legend.title = element_text(size = s1), 
        plot.title = element_text(size=s1),  axis.title.y=element_text(size=s1),
        axis.text.x = element_text(angle = 0, vjust=0.5),
        axis.text.y = element_text(angle = 90, vjust=0.5, hjust = 0.5),legend.text = element_text(size = s1), 
        legend.key = element_rect(fill=NA)) +
  geom_boxplot(color="black", fill=c(myBlue,myRed))

new3

################################################################################
## NEW NEW PLOTS
################################################################################

dfCntByRef<-aggregate(ID~ref+Type, data=df0, FUN=length1M) 
dfCntByRef$ref<-gsub("mRatBN7", "mRatBN7.2", dfCntByRef$ref)
dfCntByRef$ref<-gsub("Rnor6", "Rnor_6.0", dfCntByRef$ref)
dfCntByRef
#         ref  Type       ID
# 1 mRatBN7.2 Indel 1.615870
# 2  Rnor_6.0 Indel 3.527568
# 3 mRatBN7.2   SNP 5.088144
# 4  Rnor_6.0   SNP 8.286401


pSNP<-ggplot(data=subset(dfCntByRef, Type=="SNP"), 
              aes(x=ref, y=ID, fill=ref))+ 
              theme_bw()+
              geom_bar(stat="identity")+
              xlab("")+ylab("Count, x 1M")+
              theme(legend.position = "none", legend.title = element_text(size = s1), 
                    plot.title = element_text(size=s1),  axis.title.y=element_text(size=s1),
                    axis.text.x = element_text(angle = 0, vjust=0.5),
                    axis.text.y = element_text(angle = 90, vjust=0.5, hjust = 0.5),legend.text = element_text(size = s1), 
                    legend.key = element_rect(fill=NA))+
              ggtitle("Total number of SNPs")+
              ylim(c(0,10))

pIndel<-ggplot(data=subset(dfCntByRef, Type=="Indel"), aes(x=ref, y=ID, fill=ref))+ 
              theme_bw()+
              geom_bar(stat="identity")+
              xlab("")+ylab("Count, x 1M")+
              theme(legend.position = "none", legend.title = element_text(size = s1), 
                    plot.title = element_text(size=s1),  axis.title.y=element_text(size=s1),
                    axis.text.x = element_text(angle = 0, vjust=0.5),
                    axis.text.y = element_text(angle = 90, vjust=0.5, hjust = 0.5),legend.text = element_text(size = s1), 
                    legend.key = element_rect(fill=NA))+
              ggtitle("Total number of Indels")+
              ylim(c(0,10))

pSNP+pIndel



#ggarrange(plot1,               
#  ggarrange(p.AltSNP0, altsnpvio0, ncol = 2,nrow = 1, labels = c("D", "E")),
#  plot2,
#  ggarrange(p.AltIndel0,  altindelvio0, ncol = 2,nrow = 1, labels = c("F", "G")),
#  plot3,
#  ggarrange(p.HETSNP0, p.HETIndel0, ncol = 2,nrow = 1, labels = c("H", " I")),
#  nrow = 3,ncol = 2,
#  labels = c("A","","B","","C")       # Label of the line plot
#) 

########################################
## NEW ARRANGEMENT                    ##
########################################

## Getting somewhere, do not delete.
#ggarrange(ggarrange(new2, new3, ncol = 2,nrow = 1, labels = c("A", "B")),       
#          ggarrange(p.AltSNP0, altsnpvio0, ncol = 2,nrow = 1, labels = c("C", "D")),
#          ggarrange(p.AltIndel0,  altindelvio0, ncol = 2,nrow = 1, labels = c("E", "F")),
#          ggarrange(p.HETSNP0, p.HETIndel0, ncol = 2,nrow = 1, labels = c("G", " H")),
#          nrow = 4,ncol = 1,
#          labels = c("","","","","")       # Label of the line plot


## Horizontal version
#ggarrange(  
#  ggarrange( new3, new2, p.AltSNP0, p.AltIndel0, ncol = 4,nrow = 1, labels = c("A", "B","C","D")),
#  ggarrange( altsnpvio0, altindelvio0,  p.HETSNP0, p.HETIndel0, ncol = 4,nrow = 1, labels = c("E", "F","G","H")),
#  nrow = 2,ncol = 1,
#  labels = c("","","","","")       # Label of the line plot
#  ) 

pdf(file="figure2_new_horizontal.pdf", width=9, height=5, paper="special")

## Horizontal version
ggarrange(  
  ggarrange( new3,  pSNP,  altsnpvio0,  ncol = 3,nrow = 1, labels = c("A", "C","E")),
  ggarrange( new2, pIndel, altindelvio0, ncol = 3,nrow = 1, labels = c("B","D", "F")),
  nrow = 2,ncol = 1,
  labels = c("","","","","")       # Label of the line plot
) 

dev.off()

pdf(file="figure2_new_vertical.pdf", width=5, height=8, paper="special")

ggarrange(  
  ggarrange( new3, new2,  ncol = 2,nrow = 1, labels = c("A", "B")),
  ggarrange( pSNP, pIndel, ncol = 2,nrow = 1, labels = c("C","D")),
  ggarrange( altsnpvio0, altindelvio0, ncol = 2,nrow = 1, labels = c("E","F")),
  nrow = 3,ncol = 1,
  labels = c("","","","","")       # Label of the line plot
) 



dev.off()

pdf(file="supplementary_figure2_alt.pdf", width=8, height=4, paper="special")

ggarrange(HETsnpvio0, HETindelvio0,
          nrow = 1,ncol = 2,
          labels = c("A","B")       # Label of the line plot
) 

dev.off()

pdf(file="supplementary_figure2b_alt.pdf", width=6, height=6, paper="special")

ggarrange(p.comAltSNP0,
          p.comAltIndel0,
          p.comHETSNP0,
          p.comHETIndel0,
          nrow = 2,ncol = 2,
          labels = c("A","B","C","D")       # Label of the line plot
) 

dev.off()


ggplot(df0_both, aes(x=reorder(ref,gems_detected), y=gems_detected, fill=ref)) + geom_boxplot(fill="white") + xlab("ref")



#memory.limit()


## For the paper:

# how many RN6 SNP and indel?
table(rn6$Type[rn6$Qual >=30])

# how many BN7 SNP and indel?
table(bn7$Type[bn7$Qual >=30])









