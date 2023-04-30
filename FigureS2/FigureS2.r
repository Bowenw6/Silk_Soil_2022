
library(ggplot2)
library(ggpubr)
wbwPalette <- c("#55AFFF", "#F4B800")

#16S_α_diversity
sum_table_16S <- read.table('16S_summary.alpha_diversity.txt', sep = '\t', row.name = 1, header = TRUE)
boxdat <- as.data.frame(sum_table_16S)


compaired <- list(c('Silk', 'Soil'))

#Shannon index
  p1<- ggplot(boxdat, aes(Group, shannon)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Shannon index")+xlab(NULL)+ 
    theme(legend.position="none") 
 
  
  
  
  
#Chao1

p2<- ggplot(boxdat, aes(Group, chao)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Chao1 index")+xlab(NULL)+ 
  theme(legend.position="none") 





#Simpson index
p3<- ggplot(boxdat, aes(Group, simpson)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Simpson index")+xlab(NULL)+ 
  theme(legend.position="none") 






#Ace index
p4<- ggplot(boxdat, aes(Group, ace)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) +
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+
  scale_fill_manual(values=wbwPalette) + 
  ylab("Ace index")+xlab(NULL)+ 
  theme(legend.position="none") 





library(patchwork)
pp1<- p1+p2+p3+p4+plot_layout(ncol=5, byrow= TRUE) + 
  theme(legend.position = c(1.6,0.5), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=12), 
        legend.key.size = unit(20, "pt")) 
plot(pp1)
#save image 3inchs*10inchs









#rm(list = ls())









#ITS_α_diversity
sum_table_ITS <- read.table('ITS_summary.alpha_diversity.txt', sep = '\t', row.name = 1, header = TRUE)
boxdat1 <- as.data.frame(sum_table_ITS)


compaired <- list(c('Silk', 'Soil'))


p11<- ggplot(boxdat1, aes(Group, shannon)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Shannon index")+xlab(NULL)+ 
  theme(legend.position="none") 





#Chao1


p22<- ggplot(boxdat1, aes(Group, chao)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Chao1 index")+xlab(NULL)+ 
  theme(legend.position="none") 



#Simpson index
p33<- ggplot(boxdat1, aes(Group, simpson)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Simpson index")+xlab(NULL)+ 
  theme(legend.position="none") 




#Ace index
p44<- ggplot(boxdat1, aes(Group, ace)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) +
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) +
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Ace index")+xlab(NULL)+ 
  theme(legend.position="none") 



library(patchwork)
pp2<-p11+p22+p33+p44+plot_layout(ncol=5, byrow= TRUE) + 
  theme(legend.position = c(1.6,0.5),
        legend.text=element_text(size=10), 
        legend.title=element_text(size=12), 
        legend.key.size = unit(20, "pt")) 
plot(pp2)
#save image 3inchs*10inchs



#Annotation from ggplot2 stat_compare_means function
####In other words, we use the following convention for symbols indicating statistical significance:
#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001


