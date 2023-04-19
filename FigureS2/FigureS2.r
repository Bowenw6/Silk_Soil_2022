
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
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + #设置一下误差线
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + #离群值是黑色，与点图重合，所以此处颜色调节为NA
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + #vjust是调节星标的垂直位置
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + #去掉各种背景
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ #shape=21是空心圆点
  scale_fill_manual(values=wbwPalette) + #要用scale_fill_manual改jitter的颜色scale_color_manual无效
  ylab("Ace index")+xlab(NULL)+ #加y轴名称，隐藏x轴名称
  theme(legend.position="none") #隐藏图例





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
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + #设置一下误差线
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + #离群值是黑色，与点图重合，所以此处颜色调节为NA
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "wilcox.test", vjust = 0.6) + #vjust是调节星标的垂直位置
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + #去掉各种背景
  geom_jitter(aes(fill = Group), shape=21, width = 0.2, size = 2, alpha = 0.8)+ #shape=21是空心圆点
  scale_fill_manual(values=wbwPalette) + #要用scale_fill_manual改jitter的颜色scale_color_manual无效
  ylab("Ace index")+xlab(NULL)+ #加y轴名称，隐藏x轴名称
  theme(legend.position="none") #隐藏图例



library(patchwork)
pp2<-p11+p22+p33+p44+plot_layout(ncol=5, byrow= TRUE) + 
  theme(legend.position = c(1.6,0.5), #图例位置
        legend.text=element_text(size=10), #Silk和Soil字体大小
        legend.title=element_text(size=12), #Group字体大小
        legend.key.size = unit(20, "pt")) #整个图例的大小
plot(pp2)
#save image 3inchs*10inchs



#Annotation from ggplot2 stat_compare_means function
####In other words, we use the following convention for symbols indicating statistical significance:
#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001


