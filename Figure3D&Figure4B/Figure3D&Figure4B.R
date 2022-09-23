rm(list=ls())

library(RColorBrewer)
display.brewer.all()
brewer.pal(12,'Paired')



#Figure3D
data_16S <- read.table('all_16S.txt', sep = '\t', row.names = 1,header = T)
#Top 1582 row, which total_tags >= 50
data_16S<-data_16S[1:1582,]




bk = unique(c(seq(-1.2,1.2, length=50)))


group_siso<- read.table("group.txt", sep='\t',row.names = 1,  header=FALSE)
colnames(group_siso)<-c("Groups")

ann_colors<- list(Groups = c(Silk="#55AFFF",Soil="#F4B800"))

#plot
library(pheatmap)
p1<- pheatmap(data_16S, 
              scale = 'row',
              cluster_cols = T,
              breaks = bk,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = FALSE,
              color = colorRampPalette(c("#6A3D9A", 
                                         "white",
                                         "#E31A1C"))(50), 
              legend = T,
              border_color = 'NA',
              treeheight_row=20,
              treeheight_col=20,
              fontsize = 12.5,
              fontsize_col = 12.5,
              angle_col=45,
              annotation_col = group_siso, #注释的数据来源
              annotation_colors = ann_colors, #注释的颜色
              annotation_names_col = FALSE) #在图中不显示注释的名字



#save pdf 4*10 inches












#Figure4B
data_ITS <- read.table('all_ITS.txt', sep = '\t', row.names = 1,header = T)
#Top 742 row, which total_tags >= 50
data_ITS<-data_ITS[1:742,]



bk2 = unique(c(seq(-1.2,1.2, length=50)))


group_siso<- read.table("group.txt", sep='\t',row.names = 1,  header=FALSE)
colnames(group_siso)<-c("Groups")

ann_colors<- list(Groups = c(Silk="#55AFFF",Soil="#F4B800"))

#plot
library(pheatmap)
p2<- pheatmap(data_ITS, 
              scale = 'row',
              cluster_cols = T,
              breaks = bk2,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = FALSE,
              color = colorRampPalette(c("#6A3D9A", 
                                         "white",
                                         "#E31A1C"))(50), 
              legend = T,
              border_color = 'NA',
              treeheight_row=20,
              treeheight_col=20,
              fontsize = 12.5,
              fontsize_col = 12.5,
              angle_col=45,
              annotation_col = group_siso,
              annotation_colors = ann_colors,
              annotation_names_col = FALSE)


#save pdf 4*10 inches






