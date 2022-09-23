#BugBase
library(ggplot2)
library(ggpubr)
wbwPalette <- c("#55AFFF", "#F4B800")


BugBase_dat <- read.table('predictions.txt', sep = '\t', row.name = 1, header = TRUE)
boxdat <- as.data.frame(BugBase_dat)


compaired <- list(c('Silk', 'Soil'))


  p1<- ggplot(boxdat, aes(Group, Aerobic)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Aerobic")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 

  
  p2<- ggplot(boxdat, aes(Group, Facultatively_Anaerobic)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Facultatively Anaerobic")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  

  p3<- ggplot(boxdat, aes(Group, Anaerobic)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Anaerobic")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  

  p4<- ggplot(boxdat, aes(Group, Contains_Mobile_Elements)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "t.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Contains Mobile Elements")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  
  

  p5<- ggplot(boxdat, aes(Group, Gram_Positive)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Gram Positive")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 

  

  p6<- ggplot(boxdat, aes(Group, Gram_Negative)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Gram Negative")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 

  

  p7<- ggplot(boxdat, aes(Group, Stress_Tolerant)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.1) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Stress Tolerant")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 

  

  p8<- ggplot(boxdat, aes(Group, Forms_Biofilms)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Forms Biofilms")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 

  

  p9<- ggplot(boxdat, aes(Group, Potentially_Pathogenic	)) + 
    stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
    geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
    stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.1) + 
    theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
    geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
    scale_fill_manual(values=wbwPalette) + 
    ylab("Relative Abundance")+
    xlab(NULL)+ 
    labs(title = "Potentially Pathogenic")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") 
  
  
#Merge 9 pics
library(patchwork)
p1 + p2 + p3 + p5 + p6 + p4 + p8 + p7 + p9 + plot_layout(ncol=5, byrow= TRUE) + 
  theme(legend.position = c(1.7,0.5), 
        legend.text=element_text(size=20), 
        legend.title=element_text(size=35), 
        legend.key.size = unit(40, "pt")) 
  

#save image 8*12 inches
  

