#16S
library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)



raw_data <- read.table("16S_pathway.L2.txt", header = TRUE, sep="\t", row.names = 1)
groups <- read.table("group.txt",sep = "\t",header = F,colClasses = c("character"))
groups <- as.list(groups)

data<- as.data.frame(lapply(raw_data, function(x)x / sum(x)))

row.names(data)<- row.names(raw_data)

data <- t(data)
data[is.na(data)] <- 0
data <- vegdist(data,method = "bray")
pcoa<- pcoa(data, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","group")
pich=c(21:24)
wbwPalette <- c("#55AFFF", "#F4B800")
Palette <- c("#000000", "#000000")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

otu.adonis=adonis(data~V2,data = groups,distance = "bray") 

ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(colour=group,shape=group,fill=group),size=12)+
  geom_text(aes(x = 0.005,y = -0.035,label = paste("PERMANOVA:\n    Silk VS Soil\n    p-value < ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),size = 12,hjust = 0)+ 
  stat_ellipse(aes(fill = group),geom = "polygon",level = 0.95,alpha = 0.2)+
  scale_shape_manual(values=pich)+
  scale_colour_manual(values=Palette)+
  scale_fill_manual(values=wbwPalette)+
  labs(title="PCoA - The function of microbiome - Bacteria") + 
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  
  ylim(-0.055,0.025)+
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34),
        axis.title.y=element_text(colour='black', size=34),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_blank(),
        legend.text=element_text(size=34),
        legend.key=element_blank(),legend.position = c(0.9,0.88),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1.6,"cm"))+
  theme(plot.title = element_text(size=34,colour = "black",hjust = 0.5,face = "bold"))
#save image 12*16 inches







rm(list=ls())







#ITS
library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)

raw_data <- read.table("ITS_path_abun_unstrat_descrip.txt", header = TRUE, sep="\t", row.names = 1)
groups <- read.table("group.txt",sep = "\t",header = F,colClasses = c("character"))
groups <- as.list(groups)

data<- as.data.frame(lapply(raw_data, function(x)x / sum(x)))

row.names(data)<- row.names(raw_data)

data <- t(data)
data[is.na(data)] <- 0
data <- vegdist(data,method = "bray")
pcoa<- pcoa(data, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups$V2)
colnames(plotdata) <-c("sample","PC1","PC2","group")
pich=c(21:24)
wbwPalette <- c("#55AFFF", "#F4B800")
Palette <- c("#000000", "#000000")
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)

otu.adonis=adonis(data~V2,data = groups,distance = "bray") 

ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(colour=group,shape=group,fill=group),size=12)+
  geom_text(aes(x = 0.05,y = -0.06,label = paste("PERMANOVA:\n    Silk VS Soil\n    p-value = ",otu.adonis$aov.tab$`Pr(>F)`[1],sep = "")),size = 12,hjust = 0)+ 
  stat_ellipse(aes(fill = group),geom = "polygon",level = 0.95,alpha = 0.2)+
  scale_shape_manual(values=pich)+
  scale_colour_manual(values=Palette)+
  scale_fill_manual(values=wbwPalette)+
  labs(title="PCoA - The function of microbiome - Fungi") + 
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  
  ylim(-0.07,0.055)+ 
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34),
        axis.title.y=element_text(colour='black', size=34),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_blank(),
        legend.text=element_text(size=34),
        legend.key=element_blank(),legend.position = c(0.9,0.88),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1.6,"cm"))+
  theme(plot.title = element_text(size=34,colour = "black",hjust = 0.5,face = "bold"))
#save image 12*16 inchs

