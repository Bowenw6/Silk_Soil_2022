#生成进化树的注释文件


#筛选ITS测序的Enriched和Depleted的ASVs
#该表格来源于edgeR分析的结果
diff_result_ITS<-read.csv("diff_result_ITS.csv", row.names = 1,na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA
filtered_diff_data_ITS<- diff_result_ITS[,c("logCPM","threshold","Phylum")]
#把未分类注释到门水平的ASVs的名称改为Unknown
filtered_diff_data_ITS[is.na(filtered_diff_data_ITS)]<- "Unknown"
#选出Enriched和Depleted的ASVs
Enriched_data_ITS<- filtered_diff_data_ITS[filtered_diff_data_ITS$threshold=="Enriched on silk",]
Depleted_data_ITS<- filtered_diff_data_ITS[filtered_diff_data_ITS$threshold=="Depleted on silk",]
#拼接两个表格
Signif_diff_data_ITS<-rbind(Enriched_data_ITS,Depleted_data_ITS)
#也可以用下面的代码，直接筛选threshold不等于Not significant的(这个命令其实更好...)
#Signif_diff_data_ITS<-filtered_diff_data_ITS[filtered_diff_data_ITS$threshold!="Not significant",]


#提取差异显著的ASVs的代表性序列
Fasta_ITS<-read.table("ITS_all.otus.representative.fasta")
re_Fasta_ITS<-Fasta_ITS[seq(1,nrow(Fasta_ITS),2),]#提取奇数行放在新表格的V1
re_Fasta_ITS<-as.data.frame(re_Fasta_ITS)
re_Fasta_ITS$V2<-Fasta_ITS[seq(0,nrow(Fasta_ITS),2),]#提取偶数行放在新表格的V2
colnames(re_Fasta_ITS)<- c("ASVs_names","Sequences")#改个列名
re_Fasta_ITS$ASVs_names<- gsub(">A*","A",re_Fasta_ITS$ASVs_names) #把ASV编号前的>去掉
row.names(re_Fasta_ITS)<- re_Fasta_ITS$ASVs_names


#生成差异ASVs的fasta序列
#把re_Fasta_ITS表按照Signif_diff_data_ITS表行名进行匹配筛选
Signif_diff_Fasta_ITS<-re_Fasta_ITS[match(row.names(Signif_diff_data_ITS), row.names(re_Fasta_ITS)),]
#检验一下新生成的表格和Signif_diff_data_ITS的行名顺序是否匹配
#print(row.names(Signif_diff_Fasta_ITS)==row.names(Signif_diff_data_ITS))
#把基因名前面加上>以使之成为fasta格式
Signif_diff_Fasta_ITS$ASVs_names<- gsub("ASV",">ASV",Signif_diff_Fasta_ITS$ASVs_names)
#导出fasta文件，列名行名和双引号都不要
write.table(Signif_diff_Fasta_ITS, "Signif_diff_Fasta_ITS.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
#重新读取文件
text_Fasta_ITS<-readLines("Signif_diff_Fasta_ITS.fasta",encoding = "UTF-8")
#把逗号换为回车
text_Fasta_ITS<-gsub(",","\n",text_Fasta_ITS)
#重新导出最终fasta格式文件
writeLines(text_Fasta_ITS,"Signif_diff_Fasta_ITS.fasta")


#in MacOS command line, perform alignment with muscle
#muscle -in Signif_diff_Fasta_ITS.fasta -out aligned_Signif_diff_Fasta_ITS.fasta 
#in MEGA software, output tree file "Signif_diff_tree_ITS.nwk" (第一次出图用的这个树)
# in server, construct tree with iqtree2（第二次再出图，这个软件更合理）
#nohup iqtree -s aligned_Signif_diff_Fasta_ITS.fasta -m MFP -B 1000 --bnni -T AUTO > Signif_ITS_nohup.out & 


#把表达量总表筛选出差异ASVs用于后续热图绘制
all_ASVs_ITS<- read.table(file = "ITS_all.otus.profiling.txt",header = TRUE, row.names = 1, sep = "\t", na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA) 
Signif_diff_ASVs_ITS<-all_ASVs_ITS[match(row.names(Signif_diff_data_ITS), row.names(all_ASVs_ITS)),] #保留差异显著的ASVs的丰度数据
Signif_diff_ASVs_ITS$V32<- rownames(Signif_diff_ASVs_ITS) #把第32列改为行名，即ASVs的ID
colnames(Signif_diff_ASVs_ITS)<- gsub("V32","ASVs", colnames(Signif_diff_ASVs_ITS)) #把colnames中的V32改成ASVs
Signif_diff_ASVs_ITS[is.na(Signif_diff_ASVs_ITS)]<- "Unknown" #把未注释到门水平的命名为“Unknown”
Signif_diff_ASVs_ITS<-Signif_diff_ASVs_ITS[,c(32, 1:31)] #把ASVs放在前面美观一些
#检验一下Signif_diff_ASVs_ITS和Signif_diff_data_ITS的行名是否一致
#print(row.names(Signif_diff_ASVs_ITS)==row.names(Signif_diff_data_ITS))
Signif_diff_ASVs_ITS<-cbind(Signif_diff_ASVs_ITS,Signif_diff_data_ITS)
Signif_diff_ASVs_ITS<-Signif_diff_ASVs_ITS[,1:34] #把第35列多余的Phylum去掉
write.table(Signif_diff_ASVs_ITS,"Signif_diff_ASVs_ITS.txt", sep = "\t", row.names = FALSE, col.names =TRUE, quote =FALSE) #输出表格，保留列名不要行名，不要双引号

#Signif_diff_ASVs_ITS.txt文件即为进化树的注释文件




