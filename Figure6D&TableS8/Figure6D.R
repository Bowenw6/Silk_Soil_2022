library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(ggnewscale)
library(ggClusterNet)
library(ggstar)
library(treeio)

tree_ITS<-read.tree("aligned_Signif_diff_fasta_ITS.fasta.treefile")
data_ITS<-fortify(tree_ITS)
map_ITS<- read.table("Signif_diff_ASVs_ITS.txt", header = TRUE, sep = "\t")


p<-ggtree(tree_ITS, aes(col=Phylum), layout="circular", size=0.1) %<+% map_ITS + 
  geom_tippoint(aes(color=Phylum, shape=threshold), size=3, alpha =0.8)+
  scale_shape_manual(values=c("Enriched on silk" = 1, 
                              "Depleted on silk" = 4))+
  geom_tiplab(aes(label=NA, col=Phylum), hjust=-0.5, align=TRUE, linesize=0.2, alpha=0.5)+
  geom_tiplab(aes(label=NA, col=NA), size=0.8)+
  theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) 

p<- open_tree(p, 90) %>% rotate_tree(90) 



bar_ITS<-as.data.frame(map_ITS[,c(1,33,34)], row.names = map_ITS$ASVs)
colnames(bar_ITS)<- c("names","number_log_CPM","log_CPM")
bar_ITS$names<-as.factor(bar_ITS$names)
bar_ITS$phylum_name<-as.factor(bar_ITS$number_log_CPM)

p1<-p + 
  geom_fruit(data=bar_ITS, geom=geom_bar, mapping=aes(y=names, x=number_log_CPM, fill=c("white")), orientation="y", 
             stat="identity", axis.params = NA) + 
  scale_fill_manual(values=c("white"))+ 
  new_scale_fill()+ 
  geom_fruit(data=bar_ITS, geom=geom_bar, mapping=aes(y=names, x=number_log_CPM, fill=log_CPM), orientation="y",
             stat="identity", axis.params = list(axis="x", vjust= 0.1, text.size=3))+
  scale_fill_manual(values=c("#F4B800","#55AFFF"))+
  new_scale_fill()





tmp<- read.table("Signif_diff_ASVs_ITS.txt",header = TRUE, row.names = 1, sep = "\t")
tmp<- tmp[,1:30] 
tmp$Silk<- apply(tmp[,1:15],1,mean) 
tmp$Soil<- apply(tmp[,16:30],1,mean) 
tmp<-tmp[,31:32] 
tmp<-log10(tmp+1) 



p2<- gheatmap(p1, tmp, offset = 0, width=0.15, font.size=6, color = NULL, hjust = -0.1, colnames_level=colnames(tmp), colnames_angle=0, legend_title=" log10 (Mean abundance)",
              low="palegreen3", high ="darkorange3") + theme(legend.position = c(0.5,0.75))


p2
#save image 12inches*inches
#rearrange legends with AI



