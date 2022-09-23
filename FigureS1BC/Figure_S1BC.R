
rm(list=ls())
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(12,'Set3')
brewer.pal(12,'Set3')
Palette<-c(brewer.pal(4,'Set3'))


#FigureS1B
abundance <- read.table("EDS.txt",header = TRUE, row.names = 1, sep = "\t")
abundance <- as.matrix(abundance)

f.abundance <-abundance

CT_EDS_f.abundance<-f.abundance/100



library(reshape2)
taxon <- melt(CT_EDS_f.abundance)
colnames(taxon) <- c("Element","variable","value")
library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Element, stratum = Element))
p1 <- p + geom_alluvium(aes(fill = Element),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Element),width = 0.6)
p2 <- p1 + ylab(label = "Relative contents") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
#p4 <- p3 
 theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
 theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+
 theme(axis.title.y = element_text(size = 18,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Relative contents of elements based on EDS")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
  #theme(text = element_text(family = "Times"))
p5
#导出PDF_5*8inchs








#FigureS1C
taxon$Groups<-ifelse(substr(taxon$variable,start = 1,stop = 2)=="CK","CK","Treatment") #前两位字符是CK的Group写为CK，否则为Treatment



library(ggplot2)
library(ggpubr)
wbwPalette <- c("#55AFFF", "#F4B800")



compaired <- list(c('CK', 'Treatment'))


p6<- ggplot(taxon, aes(Groups, value)) + 
  stat_boxplot(geom = 'errorbar', width = 0.4, size=0.3) + 
  geom_boxplot(fill= "gray", width = 0.6, outlier.colour = NA, size=0.8) + 
  stat_compare_means(comparisons = compaired, label = "p.signif", method = "t.test", vjust = 0.6) + 
  theme_bw() + theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) + 
  geom_jitter(aes(fill = Groups), shape=21, width = 0.2, size = 2, alpha = 0.8)+ 
  scale_fill_manual(values=wbwPalette) + 
  ylab("Relative contents")+xlab(NULL)+ 
  theme(legend.position="none") + 
  facet_wrap("Element",scales = "free") +
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 18,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  theme(strip.text = element_text(size = 18))


p6


#Merge 2 pics
library(patchwork)
p7<-p5+p6+plot_layout(widths = c(16,9))

p7

#save images,9inchs*20inchs









