#DOC plot Reference:Fungal community assembly in drought-stressed sorghum shows stochasticity, selection, and universal ecological dynamics, Nature communication 

rm(list=ls())
library(ggplot2)
library(vegan)
library(Hmisc) 

ASV_all <- read.delim("all_16S.txt", header = TRUE ,row.names=1)


A=ASV_all 
C=A/rowSums(A)
ASV_all_1<-t(C)

ASV_all_1[ASV_all_1>0]<-1

ASV_all_1<-t(ASV_all_1)
ASV_filtered_16S<-ASV_all[which(rowSums(ASV_all_1)>15),]



summary(colSums(ASV_filtered_16S))

ASV_Silk_16S<- ASV_filtered_16S[,1:15]
ASV_Soil_16S<- ASV_filtered_16S[,16:30]


#16S_Silk group DOC plot
library(DOC)
Silk_16S_DOC_results <- DOC(ASV_Silk_16S, R = 1000)
#The P-value was calculated as the fraction of bootstrap runs (n_boot=1000) resulting with non-negative slopes.
silk_p_value<-(sum(Silk_16S_DOC_results$LME$Slope >= 0) + 1) / (length(Silk_16S_DOC_results$LME$Slope) + 1)
p1<-plot(Silk_16S_DOC_results) + theme_bw() +geom_text(aes(x = 0.45,y = 0.95,
                                          label = paste("Silk-Prokaryotes", '\n', "P=", round(silk_p_value, 4)))) 
p1
#save as pdf 4inches*4inches


#16S_Soil group DOC plot
library(DOC)
Soil_16S_DOC_results <- DOC(ASV_Soil_16S, R = 1000)
#The P-value was calculated as the fraction of bootstrap runs (n_boot=1000) resulting with non-negative slopes.
soil_p_value<-(sum(Soil_16S_DOC_results$LME$Slope >= 0) + 1) / (length(Soil_16S_DOC_results$LME$Slope) + 1)
p2<-plot(Soil_16S_DOC_results) + theme_bw() +geom_text(aes(x = 0.85,y = 0.95,
                                                       label = paste("Soil-Prokaryotes", '\n', "P=", round(soil_p_value, 4)))) #调节x和y调节字体位置
p2
#save as pdf 4inches*4inches

#Merge two pics
library(ggpubr)
p3<-ggarrange(p1, p2, ncol = 2, nrow = 1)
p3
#save as pdf 4inches*8inches







#Boxplot of beta_nti and RCbray
library(ggplot2)
library(ggpubr)
wbwPalette <- c("#55AFFF", "#F4B800")


beta_nti_Silk <- read.table("beta_nti_silk_melted.txt", sep = '\t', row.name = 1, header = TRUE)
beta_nti_Soil <- read.table("beta_nti_soil_melted.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_silk<- read.table("raup_crick_silk.txt", sep = '\t', row.name = 1, header = TRUE)
raup_crick_soil<- read.table("raup_crick_soil.txt", sep = '\t', row.name = 1, header = TRUE)


beta_nti_Silk$Group<-rep(c("Silk"),times=dim(beta_nti_Silk)[1])
beta_nti_Soil$Group<-rep(c("Soil"),times=dim(beta_nti_Soil)[1])
sum_beta_nti<- rbind(beta_nti_Silk[,3:4],beta_nti_Soil[,3:4])
colnames(sum_beta_nti)<-c("betaNTI", "Group")

raup_crick_silk$Group<-rep(c("Silk"),times=dim(raup_crick_silk)[1])
raup_crick_soil$Group<-rep(c("Soil"),times=dim(raup_crick_soil)[1])
sum_raup_crick<- rbind(raup_crick_silk[,3:4],raup_crick_soil[,3:4])
colnames(sum_raup_crick)<-c("RCbray", "Group")

compaired <- list(c('Silk', 'Soil'))
#################


#betaNTI boxplot

p4<- ggplot(sum_beta_nti, aes(Group, betaNTI)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= c("#55AFFF", "#F4B800"), width = 0.6, outlier.colour = NA, size=0.8, alpha=0.3) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 1.2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("betaNTI")+xlab(NULL)+ 
  ylim(-3,7.5)+ 
  geom_hline(aes(yintercept=2), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=-2), linetype='dashed',colour="red",alpha = 0.8) +
  theme(legend.position="none") 

p4

#save image 4inchs*3inchs


#RCbray plot

p5<- ggplot(sum_raup_crick, aes(Group, RCbray)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= c("#55AFFF", "#F4B800"), width = 0.6, outlier.colour = NA, size=0.8, alpha=0.3) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = -0.2) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 1.2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("RCbray")+xlab(NULL)+ 
  ylim(-2.1,2.1)+ 
  geom_hline(aes(yintercept=0.95), linetype='dashed',colour="red",alpha = 0.8) +
  geom_hline(aes(yintercept=-0.95), linetype='dashed',colour="red",alpha = 0.8) +
  theme(legend.position="none") 

p5
#save image 4inchs*3inchs


#Merge two pics 
library(ggpubr)
p6<-ggarrange(p4, p5, ncol = 2, nrow = 1)
p6 

#save pdf,4inches*6inches


