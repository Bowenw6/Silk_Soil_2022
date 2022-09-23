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

#开始画图
p<-ggtree(tree_ITS, aes(col=Phylum), layout="circular", size=0.1) %<+% map_ITS + # %<+% map # 引入注释文件
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tippoint(aes(color=Phylum, shape=threshold), size=3, alpha =0.8)+
  scale_shape_manual(values=c("Enriched on silk" = 1, 
                              "Depleted on silk" = 4))+
  #加树和外圈之间的连线
  geom_tiplab(aes(label=NA, col=Phylum), hjust=-0.5, align=TRUE, linesize=0.2, alpha=0.5)+
  # 端点颜色、大小
  geom_tiplab(aes(label=NA, col=NA), size=0.8)+
  # 注释、注释的颜色
  theme(legend.title=element_text(face="bold"), legend.position="right", legend.box="horizontal", legend.text=element_text(size=rel(0.7))) 
  #外圈添加柱状图
  # 图例位置、文字大小
  #xlim(NA, max(data_ITS$x)*1.3)

p<- open_tree(p, 90) %>% rotate_tree(90) #树开口90度，再把树旋转90度


#外圈绘制柱状图
bar_ITS<-as.data.frame(map_ITS[,c(1,33,34)], row.names = map_ITS$ASVs)
colnames(bar_ITS)<- c("names","number_log_CPM","log_CPM")
bar_ITS$names<-as.factor(bar_ITS$names)
bar_ITS$phylum_name<-as.factor(bar_ITS$number_log_CPM)

p1<-p + 
  geom_fruit(data=bar_ITS, geom=geom_bar, mapping=aes(y=names, x=number_log_CPM, fill=c("white")), orientation="y", #不fill
             stat="identity", axis.params = NA) + #坐标轴隐藏
  scale_fill_manual(values=c("white"))+ #颜色调为白色
  new_scale_fill()+ #以上只是加一层空白，给热图留位置
  geom_fruit(data=bar_ITS, geom=geom_bar, mapping=aes(y=names, x=number_log_CPM, fill=log_CPM), orientation="y",
             stat="identity", axis.params = list(axis="x", vjust= 0.1, text.size=3))+
  scale_fill_manual(values=c("#F4B800","#55AFFF"))+
  new_scale_fill()




#外圈加热图
#配制热图文件tmp
tmp<- read.table("Signif_diff_ASVs_ITS.txt",header = TRUE, row.names = 1, sep = "\t")
tmp<- tmp[,1:30] #只保留数值型的数据
tmp$Silk<- apply(tmp[,1:15],1,mean) #求Silk组的平均丰度
tmp$Soil<- apply(tmp[,16:30],1,mean) #求Soil组的平均丰度
tmp<-tmp[,31:32] #只保留算得的平均值
tmp<-log10(tmp+1) #将相对丰度转换为log10，+1是为了避免丰度为0无法取log

#还有另一种思路，算相对丰度之后再求log，太复杂了，需要转换两次，故不采纳了
#tmp<- t(t(tmp)/colSums(tmp,na=T))*100 #转化成相对丰度
#tmp<-as.data.frame(tmp,colnames=TRUE)
#tmp$Silk<- apply(tmp[,1:15],1,sum) #求Silk组的相对丰度
#tmp$Soil<- apply(tmp[,16:30],1,sum) #求Soil组的相对丰度
#tmp<-tmp[,31:32] #只保留算得的值
#tmp<-log10(tmp+1) #将相对丰度转换为log10

p2<- gheatmap(p1, tmp, offset = 0, width=0.15, font.size=6, color = NULL, hjust = -0.1, colnames_level=colnames(tmp), colnames_angle=0, legend_title=" log10 (Mean abundance)",
              low="palegreen3", high ="darkorange3") + theme(legend.position = c(0.5,0.75))


p2
#save image 12inches*inches
#rearrange legends with AI



