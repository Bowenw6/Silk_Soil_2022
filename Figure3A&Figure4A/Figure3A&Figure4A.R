#Figure3A
Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")
abundance <- read.table("16S_profiling.Phylum.top10.txt",header = TRUE,sep = "\t")
f.abundance <- as.data.frame(abundance)
library(reshape2)
taxon <- melt(f.abundance)
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Taxon),width = 0.6)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
#p4 <- p3 
 theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
 theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
 theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 10 phyla of bacteria")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
  #theme(text = element_text(family = "Times"))
p5
#PDF_6*18inchs





rm(list=ls())





#Figure4A
Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")
abundance <- read.table("ITS_profiling.Phylum.top10.txt",header = TRUE,sep = "\t")
f.abundance <- as.data.frame(abundance)
library(reshape2)
taxon <- melt(f.abundance)
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
library(ggalluvial)

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + 
  geom_stratum(aes(fill = Taxon),width = 0.6)
p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p3 <- p2 + scale_fill_manual(values = Palette) +
  theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm")) +
  #p4 <- p3 
  theme(axis.text.x=element_text(colour="black",size=12,face = "bold",angle = -45,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 18)) + 
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 24,face = "bold",margin = unit(c(0,1,0,1),"lines")))+
  scale_y_continuous(expand = c(0,0))
p5 <- p3 + theme(legend.text = element_text(colour = "black",size = 15)) + 
  theme(legend.title = element_text(size = 18,colour = "black")) + 
  labs(title = "Top 10 phyla of fungi")+
  theme(text=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5)) 
#theme(text = element_text(family = "Times"))
p5
#PDF_6*18inchs
