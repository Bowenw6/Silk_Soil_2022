rm(list=ls())

library(ggplot2)
library(tidyverse)
library(randomForest)
library(rfUtilities)
library(rfPermute)


All_ITS <- read.table(file = "all_ITS.txt",header = T, row.names = 1, sep = "\t")


ASV_filtered_ITS<-All_ITS[which(rowSums(All_ITS)>=30),]





ASV_filtered_ITS<-t(ASV_filtered_ITS)
ASV_filtered_ITS<-as.data.frame(ASV_filtered_ITS)

ASV_filtered_ITS$Group<-if_else(grepl("*Si", rownames(ASV_filtered_ITS)), "Silk","Soil")

ASV_filtered_ITS$Group<-as.factor(ASV_filtered_ITS$Group)



set.seed(666)

treat_rf <- randomForest(as.factor(Group) ~ ., data= ASV_filtered_ITS,importance=TRUE,proximity=TRUE)
treat_rf
plot(treat_rf, main = "RandomForest origin")
#save image 5inches*4inches

set.seed(666)
treat_perm <- rf.significance(treat_rf, ASV_filtered_ITS, nperm=99, ntree=500)
treat_perm



Group_rfP<- rfPermute(Group ~ ., data = ASV_filtered_ITS, ntree = 500,
                         na.action = na.omit, nrep = 100,num.cores = 1)
Group_dat <- importance(Group_rfP, sort.by = NULL, decreasing = TRUE)
Group_dat[,c("MeanDecreaseAccuracy","MeanDecreaseAccuracy.pval")] %>% 
  as_tibble(rownames = "names") %>% 
  mutate(label = if_else(MeanDecreaseAccuracy.pval<0.001,"***",
                         if_else(MeanDecreaseAccuracy.pval<0.01,"**", 
                                 if_else(MeanDecreaseAccuracy.pval<0.05,"*","ns")))) %>% 
  arrange(MeanDecreaseAccuracy) %>% 
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>% 
  ggplot(aes(x = names, y = MeanDecreaseAccuracy))+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = MeanDecreaseAccuracy+0.5,label = label))+
  labs(x = "", y = "% Mean decrease accuracy")+
  coord_flip()

#save image 24inches *24 inches


Group_dat<-as.data.frame(Group_dat)
Top_features<- Group_dat[which(Group_dat$MeanDecreaseAccuracy.pval<0.01),]
Top_features<- Top_features[which(Top_features$MeanDecreaseAccuracy>=2),]


write.csv(Top_features, file="Top_features_ITS.csv", quote = FALSE)

ASV_ITS_Taxonomy<-read.table(file = "ASV_ITS_Taxonomy.txt", sep = "\t", row.names = 1, header = TRUE)

ASV_ITS_Taxonomy_filtered<-ASV_ITS_Taxonomy[match(row.names(Top_features), row.names(ASV_ITS_Taxonomy)),]

ASV_ITS_Taxonomy_filtered[ASV_ITS_Taxonomy_filtered==""]<- "Unknown"

print(row.names(ASV_ITS_Taxonomy_filtered)==row.names(Top_features))

Top_features_plot<- cbind(Top_features, ASV_ITS_Taxonomy_filtered)

Top_features_plot$ASVid<-rownames(Top_features_plot)

write_csv(Top_features_plot, file = "Top_features.csv", quote = NULL)


#plot
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9,'Set3')
brewer.pal(9,'Set3')
plot_pallette<- c("#8DD3C7","#D9D9D9")

library(ggplot2)
p3<-ggplot(data= Top_features_plot, 
           aes(x=reorder(ASVid, MeanDecreaseAccuracy), y=MeanDecreaseAccuracy)) +
  geom_bar(aes(fill = Phylum), stat = "identity")+
  scale_fill_manual(values = c(plot_pallette))+
  geom_text(aes(x = reorder(ASVid, MeanDecreaseAccuracy), label = Class), hjust=1.3, size =2.8 )+ #在柱状图上加纲水平分类信息
  coord_flip() + 
  theme_classic() +
  xlab("")

p3
#save image 5*inches*8inches





#Top_features heatmap plot



All_ITS <- read.table(file = "all_ITS.txt",header = T, row.names = 1, sep = "\t")

Top_features_abundance<-All_ITS[match(row.names(Top_features_plot), row.names(All_ITS)),]

print(row.names(ASV_ITS_Taxonomy_filtered)==row.names(Top_features_plot))

ann_phylum<-Top_features_plot[,c("Phylum","ASVid")]
ann_phylum<-as.data.frame(ann_phylum[,1])
colnames(ann_phylum)<-c("Phylum")
rownames(ann_phylum)<-c(rownames(Top_features_plot))

bk = unique(c(seq(-2.5,2.5, length=50)))


group_siso<- read.table("group.txt", sep='\t',row.names = 1,  header=FALSE)
colnames(group_siso)<-c("Groups")

phylum_names<-sort(unique(Top_features_plot[,"Phylum"]))
phylum_ann_color<-plot_pallette
names(phylum_ann_color)<-phylum_names

ann_colors<- list(Groups = c(Silk="#55AFFF",Soil="#F4B800"), Phylum = phylum_ann_color)

#plot
library(pheatmap)
p4<- pheatmap(Top_features_abundance, 
              scale = 'row',
              cluster_cols = T,
              breaks = bk,
              cluster_row = T,
              show_colnames     = T,
              show_rownames     = T,
              color = colorRampPalette(c("#6A3D9A", 
                                         "white",
                                         "#E31A1C"))(50), 
              legend = T,
              border_color = 'NA',
              treeheight_row=40,
              treeheight_col=40,
              fontsize = 12.5,
              fontsize_col = 12.5,
              angle_col=45,
              annotation_col = group_siso, 
              annotation_row = ann_phylum, 
              annotation_colors = ann_colors, 
              annotation_names_col = FALSE, 
              annotation_names_row = FALSE) 
p4
#save pdf 8*12 inches


#Top_features pie chart plot  
library(dplyr)
Summarise_Top_features<-Top_features_plot %>% group_by(Phylum) %>% summarise(n = n())
library(ggplot2)
# Compute the position of labels
pie_data <- Summarise_Top_features %>% 
  arrange(desc(Phylum)) %>%  
  mutate(prop = n / sum(Summarise_Top_features$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


#plot_pie_chart

p5<- ggplot(pie_data, aes(x="", y=prop, fill=Phylum)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = n, x=1.2), color = "black", size=6) +
  scale_fill_manual(values = plot_pallette)

p5

#save pdf 3*3 inches




