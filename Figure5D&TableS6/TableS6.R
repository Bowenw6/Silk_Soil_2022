
diff_result_16S<-read.csv("diff_result_16S.csv", row.names = 1,na.strings = "") 
filtered_diff_data_16S<- diff_result_16S[,c("logCPM","threshold","Phylum")]

filtered_diff_data_16S[is.na(filtered_diff_data_16S)]<- "Unknown"

Enriched_data_16S<- filtered_diff_data_16S[filtered_diff_data_16S$threshold=="Enriched on silk",]
Depleted_data_16S<- filtered_diff_data_16S[filtered_diff_data_16S$threshold=="Depleted on silk",]

Signif_diff_data_16S<-rbind(Enriched_data_16S,Depleted_data_16S)




Fasta_16S<-read.table("16S_all.otus.representative.fasta")
re_Fasta_16S<-Fasta_16S[seq(1,nrow(Fasta_16S),2),]
re_Fasta_16S<-as.data.frame(re_Fasta_16S)
re_Fasta_16S$V2<-Fasta_16S[seq(0,nrow(Fasta_16S),2),]
colnames(re_Fasta_16S)<- c("ASVs_names","Sequences")
re_Fasta_16S$ASVs_names<- gsub(">A*","A",re_Fasta_16S$ASVs_names) 
row.names(re_Fasta_16S)<- re_Fasta_16S$ASVs_names



Signif_diff_fasta_16S<-re_Fasta_16S[match(row.names(Signif_diff_data_16S), row.names(re_Fasta_16S)),]
Signif_diff_fasta_16S$ASVs_names<- gsub("ASV",">ASV",Signif_diff_fasta_16S$ASVs_names)

write.table(Signif_diff_fasta_16S, "Signif_diff_fasta_16S.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)

text_fasta_16S<-readLines("Signif_diff_fasta_16S.fasta",encoding = "UTF-8")

text_fasta_16S<-gsub(",","\n",text_fasta_16S)

writeLines(text_fasta_16S,"Signif_diff_fasta_16S.fasta")


#in MacOS command line, perform alignment with muscle
#muscle -in Signif_diff_fasta_16S.fasta -out aligned_Signif_diff_fasta_16S.fasta 
# in server, construct tree with iqtree2
#nohup iqtree -s aligned_Signif_diff_fasta_16S.fasta -m MFP -B 1000 --bnni -T AUTO > Signif_16S_nohup.out & 



all_ASVs_16S<- read.table(file = "16S_all.otus.profiling.txt",header = TRUE, row.names = 1, sep = "\t", na.strings = "")  
Signif_diff_ASVs_16S<-all_ASVs_16S[match(row.names(Signif_diff_data_16S), row.names(all_ASVs_16S)),] 
Signif_diff_ASVs_16S$V32<- rownames(Signif_diff_ASVs_16S) 
colnames(Signif_diff_ASVs_16S)<- gsub("V32","ASVs", colnames(Signif_diff_ASVs_16S)) 
Signif_diff_ASVs_16S[is.na(Signif_diff_ASVs_16S)]<- "Unknown" 
Signif_diff_ASVs_16S<-Signif_diff_ASVs_16S[,c(32, 1:31)] 

Signif_diff_ASVs_16S<-cbind(Signif_diff_ASVs_16S,Signif_diff_data_16S)
Signif_diff_ASVs_16S<-Signif_diff_ASVs_16S[,1:34] 
#Table_S6
write.table(Signif_diff_ASVs_16S,"Signif_diff_ASVs_16S.txt", sep = "\t", row.names = FALSE, col.names =TRUE, quote =FALSE) 






