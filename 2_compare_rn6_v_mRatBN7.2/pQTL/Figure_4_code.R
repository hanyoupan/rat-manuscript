#rm(list=ls())
#setwd("D:/data/02110221_pQTL_analysis/QTL_rn6_rn7/final_code/rat_genome")
#####################################################################
# Figure 4A
#####################################################################

library(readr)
library(qvalue)
library(tidyverse)
library(plyr)
library(readxl)
###################################################################
exp <- read_excel("Data/Data_final_v2.xlsx",  sheet = "Figure4A")
###################################################################

ggplot(data = exp,  aes(x = abs(Distance),  fill = Group)) +
  geom_histogram(aes(x = abs(Distance), y = after_stat(density), alpha=0.5), 
                 position="identity", bins = 24) +
    theme_classic()

ggplot(data = exp, aes(x = abs(Distance),  col = Group)) +
  geom_density( adjust =2) +
   theme_classic()

ggplot(data = exp) +
  geom_histogram(aes(x = abs(Distance),  fill = Group, colour = Group), alpha = 0.3, position = "identity", bins = 24) +
  theme_classic()




#####################################################################
# Figure 4C
#####################################################################
library(ggplot2)
library(readxl)
exp <- read_excel("Data/Data_final_v2.xlsx",  sheet = "Figure4C")

# Calculating correlation coefficient
cor.test(exp$rn6, exp$rn7)
ggplot(data =  exp, aes( x= rn6 , y =  rn7)) +
  geom_point() + theme_classic() +
  geom_smooth(method = "lm")



