library(ggplot2)
library(edgeR)
library(ggpubr)
otu<-read.csv("ITS_ASV.csv",row.names = 1)
design<-read.table("design.txt",header = TRUE)

edgeR_enrich <- DGEList(counts=otu, 
                        group=design$Group)

edgeR_enrich <- calcNormFactors(edgeR_enrich)
otu_norm_enrich <- cpm(edgeR_enrich, normalized.lib.sizes=T, log=F)

model_mat_enrich <- model.matrix(~Group, data=design)

dge_enrich <- estimateGLMRobustDisp(edgeR_enrich, design=model_mat_enrich)

fit_enrich <- glmFit(dge_enrich, design=model_mat_enrich) 
lrt_enrich <- glmLRT(fit_enrich, coef=2) 
tt_enrich<- topTags(lrt_enrich, n=Inf, p.value=1) 
write.csv(topTags(lrt_enrich, n=Inf, p.value=1),'ITS_Silk_VS_Soil_glmLRT.csv', quote = FALSE)


tt_enrich<- as.data.frame(tt_enrich)
plotdata<- data.frame(2^tt_enrich$logCPM,
                      tt_enrich$logFC,
                      ifelse(tt_enrich$logFC>=0&tt_enrich$FDR<=0.05,"Enriched on silk",
                             ifelse(tt_enrich$logFC<0&tt_enrich$FDR<=0.05,
                                    "Depleted on silk","Not significant")))
colnames(plotdata) <- c("CPM", "logFC", "threshold")


p1<- ggplot(plotdata,aes(x=logFC,y=CPM,color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#F4B800","#55AFFF","gray"))+
  theme_bw()+
  theme(axis.line.x=element_line(linetype=1,color="black",size=0.5))+ 
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.5))+ 
 theme(legend.title = element_blank(), panel.border = element_blank())+ 
  ylab('Log2 (count per million)')+
  xlab('Log2 (fold change)')+
  coord_cartesian(ylim = c(0,15000)) 


p2<- ggplot(plotdata,aes(x=logFC,y=CPM,color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#F4B800","#55AFFF","gray"))+
  theme_bw()+
  theme(axis.line.y=element_line(linetype=1,color="black",size=0.5))+ 
  theme(legend.title = element_blank(),panel.border = element_blank())+ 
  ylab('')+
  xlab('')+
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(ylim = c(16000,32000)) +  
  scale_y_continuous(breaks = c(16000,32000,10000)) +
theme(legend.position="none") 

ggarrange(p2,p1,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")

#save image 4 inches * 5 inches

#Donat plot
Total_enrich<- tt_enrich
Total_enrich$Signif<- ifelse(Total_enrich$logFC>=0&Total_enrich$FDR<=0.05,"Enriched on silk",
                          ifelse(Total_enrich$logFC<0&Total_enrich$FDR<=0.05,
                                 "Depleted on silk","Not significant"))

Silk_enrich<- Total_enrich[Total_enrich$Signif=="Enriched on silk",]
Soil_enrich<- Total_enrich[Total_enrich$Signif=="Depleted on silk",]
Not_Signif_enrich<- Total_enrich[Total_enrich$Signif=="Not significant",]


# Create data
Donut_data <- data.frame(
  category=c("Enriched on silk", "Depleted on silk", "Not significant"),
  count=c(nrow(Silk_enrich), nrow(Soil_enrich), nrow(Not_Signif_enrich)))

# Compute percentages
Donut_data$fraction <- Donut_data$count / sum(Donut_data$count)

# Compute the cumulative percentages (top of each rectangle)
Donut_data$ymax <- cumsum(Donut_data$fraction)

# Compute the bottom of each rectangle
Donut_data$ymin <- c(0, head(Donut_data$ymax, n=-1))

# Compute label position
Donut_data$labelPosition <- (Donut_data$ymax + Donut_data$ymin) / 2

# Compute a good label
Donut_data$label <- paste0(Donut_data$category, "\n ASVs: ", Donut_data$count)


# Make the plot
ggplot(Donut_data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=6) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("#F4B800","#55AFFF","gray")) +
  scale_color_manual(values=c("#F4B800","#55AFFF","gray")) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI







Taxonomy_ITS<-read.table("ITS_Taxonomy.txt", sep="\t", header = TRUE, row.names = 1)

reorderd_tt_enrich<- tt_enrich[match(row.names(Taxonomy_ITS),row.names(tt_enrich)),]

diff_result_ITS<- cbind(reorderd_tt_enrich,Taxonomy_ITS)
diff_result_ITS$threshold<- ifelse(diff_result_ITS$logFC>=0&diff_result_ITS$FDR<=0.05,"Enriched on silk",
                                   ifelse(diff_result_ITS$logFC<0&diff_result_ITS$FDR<=0.05,
                                          "Depleted on silk","Not significant"))
write.csv(diff_result_ITS,'diff_result_ITS.csv', quote = FALSE)#把主要结果输出为文件



#Radar plot
diff_result_ITS<-read.csv("diff_result_ITS.csv", row.names = 1,na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA
radar_data<- diff_result_ITS[,c("logCPM","threshold","Phylum")]
radar_data$logCPM<- 2^radar_data$logCPM
names(radar_data)[1]<-"CPM"

Enriched_radar_data<- radar_data[radar_data$threshold=="Enriched on silk",]
Depleted_radar_data<- radar_data[radar_data$threshold=="Depleted on silk",]

Enriched_radar_data[is.na(Enriched_radar_data)]<- "Unknown"
Depleted_radar_data[is.na(Depleted_radar_data)]<- "Unknown"

Enriched_radar_data_long<-Enriched_radar_data[,c(1,3)]
Depleted_radar_data_long<-Depleted_radar_data[,c(1,3)]

Enriched_radar_plot<-as.data.frame(tapply(Enriched_radar_data_long$CPM, Enriched_radar_data_long$Phylum, sum))
Depleted_radar_plot<-as.data.frame (tapply(Depleted_radar_data_long$CPM, Depleted_radar_data_long$Phylum, sum))

Enriched_radar_plot[,2]<-Enriched_radar_plot[,1]
Enriched_radar_plot[,1]<-rownames(Enriched_radar_plot)
colnames(Enriched_radar_plot)<-c("Phylum", "CPM")
Depleted_radar_plot[,2]<-Depleted_radar_plot[,1]
Depleted_radar_plot[,1]<-rownames(Depleted_radar_plot)
colnames(Depleted_radar_plot)<-c("Phylum", "CPM")

Radar_plot<-merge(Enriched_radar_plot,Depleted_radar_plot, by="Phylum",all.x = T, all.y = T)
colnames(Radar_plot)<-c("Phylum","Enriched on silk", "Depleted on silk")
Radar_plot[is.na(Radar_plot)]<-1 #把无表达量的换成1，下一步取log2之后就变成了0


Radar_plot$`Enriched on silk`<-log2(Radar_plot$`Enriched on silk`)
Radar_plot$`Depleted on silk`<-log2(Radar_plot$`Depleted on silk`)


#Radar plot
Radar_plot$Max<- rep(17,7)  
Radar_plot$Min<- rep(0,7)   
rownames(Radar_plot)<-Radar_plot[,1] 
Radar_plot<-Radar_plot[,2:5]
Radar_plot<-Radar_plot[,c(3,4,1,2)]
library(dplyr)
Radar_plot<-Radar_plot[order(-Radar_plot$`Enriched on silk`),] #按Enriched on silk倒序排序为了画图美观
Radar_plot<-t(Radar_plot)
Radar_plot<-as.data.frame(Radar_plot)


library(fmsb)

radarchart(Radar_plot, axistype=1 , 
            #custom polygon
            pcol = c("#55AFFF", "#F4B800") , pfcol = scales::alpha(c("#55AFFF", "#F4B800"),0.3) , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            #custom labels
            vlcex=0.8)

# Add a legend
legend(x=-0.2, y=-0.2, legend = rownames(Radar_plot[-c(2,1),]), bty = "n", pch=20 , col=c("#55AFFF", "#F4B800") , text.col = "black", cex=1, pt.cex=3) #pch → marker shape

#save image as pdf 7inches*7inches
#rearrange the position of legend with AI






#Manhattan plot
diff_result_ITS<-read.csv("diff_result_ITS.csv", row.names = 1,na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA
Manhattan_data<- diff_result_ITS[,c("logCPM","FDR","threshold","Phylum")]

Manhattan_data[is.na(Manhattan_data)]<- "Unknown"

Manhattan_data$neglog10FDR = -log10(Manhattan_data$FDR) 

Manhattan_data[,6]<- rownames(Manhattan_data)
colnames(Manhattan_data)<-c("logCPM","FDR","threshold","Phylum","neglog10FDR","ASVs")

Man_top_phylum=c("Ascomycota","Mortierellomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Rozellomycota","Unknown" )
Manhattan_data[!(Manhattan_data$Phylum %in% Man_top_phylum),]$Phylum = "Low Abundance"


Manhattan_data<-Manhattan_data[order(Manhattan_data$`Phylum`),] #按Phylum排序为了画图美观

FDR_line <- min(Manhattan_data$neglog10FDR[Manhattan_data$threshold=="Depleted on silk"]) #FDR阈值线
library(ggplot2)

basic_theme <- theme(panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.line = element_line(size = 1),
                     legend.background = element_blank(),
                     legend.key = element_blank())
ggplot(Manhattan_data, aes(x=Phylum, y=neglog10FDR, color=Phylum, size=logCPM, shape=threshold)) +
  geom_hline(yintercept=FDR_line, linetype=2, color="red") +
  geom_jitter(alpha=0.7) + 
  scale_shape_manual(values=c("Enriched on silk" = 17, 
                              "Depleted on silk" = 25, 
                              "Not significant" = 20))+
  scale_size(breaks=c(6, 9, 12, 15)) +
  labs(x="ASVs", y="-log10(FDR)") +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),legend.position="right") +
  basic_theme

#save image as pdf 6 inches*12 inches

