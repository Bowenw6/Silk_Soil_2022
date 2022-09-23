#生成进化树的注释文件


#筛选16S测序的Enriched和Depleted的ASVs
#该表格来源于edgeR分析的结果
diff_result_16S<-read.csv("diff_result_16S.csv", row.names = 1,na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA
filtered_diff_data_16S<- diff_result_16S[,c("logCPM","threshold","Phylum")]
#把未分类注释到门水平的ASVs的名称改为Unknown
filtered_diff_data_16S[is.na(filtered_diff_data_16S)]<- "Unknown"
#选出Enriched和Depleted的ASVs
Enriched_data_16S<- filtered_diff_data_16S[filtered_diff_data_16S$threshold=="Enriched on silk",]
Depleted_data_16S<- filtered_diff_data_16S[filtered_diff_data_16S$threshold=="Depleted on silk",]
#拼接两个表格
Signif_diff_data_16S<-rbind(Enriched_data_16S,Depleted_data_16S)
#也可以用下面的代码，直接筛选threshold不等于Not significant的(这个命令其实更好...)
#Signif_diff_data_16S<-filtered_diff_data_16S[filtered_diff_data_16S$threshold!="Not significant",]


#提取差异显著的ASVs的代表性序列
Fasta_16S<-read.table("16S_all.otus.representative.fasta")
re_Fasta_16S<-Fasta_16S[seq(1,nrow(Fasta_16S),2),]#提取奇数行放在新表格的V1
re_Fasta_16S<-as.data.frame(re_Fasta_16S)
re_Fasta_16S$V2<-Fasta_16S[seq(0,nrow(Fasta_16S),2),]#提取偶数行放在新表格的V2
colnames(re_Fasta_16S)<- c("ASVs_names","Sequences")#改个列名
re_Fasta_16S$ASVs_names<- gsub(">A*","A",re_Fasta_16S$ASVs_names) #把ASV编号前的>去掉
row.names(re_Fasta_16S)<- re_Fasta_16S$ASVs_names


#生成差异ASVs的fasta序列
#把re_Fasta_16S表按照Signif_diff_data_16S表行名进行匹配筛选
Signif_diff_fasta_16S<-re_Fasta_16S[match(row.names(Signif_diff_data_16S), row.names(re_Fasta_16S)),]
#检验一下新生成的表格和Signif_diff_data_16S的行名顺序是否匹配
#print(row.names(Signif_diff_fasta_16S)==row.names(Signif_diff_data_16S))
#把基因名前面加上>以使之成为fasta格式
Signif_diff_fasta_16S$ASVs_names<- gsub("ASV",">ASV",Signif_diff_fasta_16S$ASVs_names)
#导出fasta文件，列名行名和双引号都不要
write.table(Signif_diff_fasta_16S, "Signif_diff_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
#重新读取文件
text_fasta_16S<-readLines("Signif_diff_fasta_16S.fasta",encoding = "UTF-8")
#把逗号换为回车
text_fasta_16S<-gsub(",","\n",text_fasta_16S)
#重新导出最终fasta格式文件
writeLines(text_fasta_16S,"Signif_diff_fasta_16S.fasta")


#in MacOS command line, perform alignment with muscle
#muscle -in Signif_diff_fasta_16S.fasta -out aligned_Signif_diff_fasta_16S.fasta 
#in MEGA software, output tree file "Signif_diff_tree_16S.nwk" (第一次出图用的这个树)
# in server, construct tree with iqtree2（第二次再出图，这个软件更合理）
#nohup iqtree -s aligned_Signif_diff_fasta_16S.fasta -m MFP -B 1000 --bnni -T AUTO > Signif_16S_nohup.out & 


#把表达量总表筛选出差异ASVs用于后续热图绘制
all_ASVs_16S<- read.table(file = "16S_all.otus.profiling.txt",header = TRUE, row.names = 1, sep = "\t", na.strings = "") #此处加na.strings = ""是为了把表格中空值替换为NA) 
Signif_diff_ASVs_16S<-all_ASVs_16S[match(row.names(Signif_diff_data_16S), row.names(all_ASVs_16S)),] #保留差异显著的ASVs的丰度数据
Signif_diff_ASVs_16S$V32<- rownames(Signif_diff_ASVs_16S) #把第32列改为行名，即ASVs的ID
colnames(Signif_diff_ASVs_16S)<- gsub("V32","ASVs", colnames(Signif_diff_ASVs_16S)) #把colnames中的V32改成ASVs
Signif_diff_ASVs_16S[is.na(Signif_diff_ASVs_16S)]<- "Unknown" #把未注释到门水平的命名为“Unknown”
Signif_diff_ASVs_16S<-Signif_diff_ASVs_16S[,c(32, 1:31)] #把ASVs放在前面美观一些
#检验一下Signif_diff_ASVs_16S和Signif_diff_data_16S的行名是否一致
#print(row.names(Signif_diff_ASVs_16S)==row.names(Signif_diff_data_16S))
Signif_diff_ASVs_16S<-cbind(Signif_diff_ASVs_16S,Signif_diff_data_16S)
Signif_diff_ASVs_16S<-Signif_diff_ASVs_16S[,1:34] #把第35列多余的Phylum去掉
#Table_S6
write.table(Signif_diff_ASVs_16S,"Signif_diff_ASVs_16S.txt", sep = "\t", row.names = FALSE, col.names =TRUE, quote =FALSE) #输出表格，保留列名不要行名，不要双引号






