
diff_result_ITS<-read.csv("diff_result_ITS.csv", row.names = 1,na.strings = "") 
filtered_diff_data_ITS<- diff_result_ITS[,c("logCPM","threshold","Phylum")]

filtered_diff_data_ITS[is.na(filtered_diff_data_ITS)]<- "Unknown"

Enriched_data_ITS<- filtered_diff_data_ITS[filtered_diff_data_ITS$threshold=="Enriched on silk",]
Depleted_data_ITS<- filtered_diff_data_ITS[filtered_diff_data_ITS$threshold=="Depleted on silk",]

Signif_diff_data_ITS<-rbind(Enriched_data_ITS,Depleted_data_ITS)


Fasta_ITS<-read.table("ITS_all.otus.representative.fasta")
re_Fasta_ITS<-Fasta_ITS[seq(1,nrow(Fasta_ITS),2),]
re_Fasta_ITS<-as.data.frame(re_Fasta_ITS)
re_Fasta_ITS$V2<-Fasta_ITS[seq(0,nrow(Fasta_ITS),2),]
colnames(re_Fasta_ITS)<- c("ASVs_names","Sequences")
re_Fasta_ITS$ASVs_names<- gsub(">A*","A",re_Fasta_ITS$ASVs_names) 
row.names(re_Fasta_ITS)<- re_Fasta_ITS$ASVs_names
Signif_diff_Fasta_ITS<-re_Fasta_ITS[match(row.names(Signif_diff_data_ITS), row.names(re_Fasta_ITS)),]
Signif_diff_Fasta_ITS$ASVs_names<- gsub("ASV",">ASV",Signif_diff_Fasta_ITS$ASVs_names)
write.table(Signif_diff_Fasta_ITS, "Signif_diff_Fasta_ITS.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)
text_Fasta_ITS<-readLines("Signif_diff_Fasta_ITS.fasta",encoding = "UTF-8")
text_Fasta_ITS<-gsub(",","\n",text_Fasta_ITS)
writeLines(text_Fasta_ITS,"Signif_diff_Fasta_ITS.fasta")


#in MacOS command line, perform alignment with muscle
#muscle -in Signif_diff_Fasta_ITS.fasta -out aligned_Signif_diff_Fasta_ITS.fasta 
# in server, construct tree with iqtree2
#nohup iqtree -s aligned_Signif_diff_Fasta_ITS.fasta -m MFP -B 1000 --bnni -T AUTO > Signif_ITS_nohup.out & 


all_ASVs_ITS<- read.table(file = "ITS_all.otus.profiling.txt",header = TRUE, row.names = 1, sep = "\t", na.strings = "")
Signif_diff_ASVs_ITS<-all_ASVs_ITS[match(row.names(Signif_diff_data_ITS), row.names(all_ASVs_ITS)),]
Signif_diff_ASVs_ITS$V32<- rownames(Signif_diff_ASVs_ITS) 
colnames(Signif_diff_ASVs_ITS)<- gsub("V32","ASVs", colnames(Signif_diff_ASVs_ITS)) 
Signif_diff_ASVs_ITS[is.na(Signif_diff_ASVs_ITS)]<- "Unknown" 
Signif_diff_ASVs_ITS<-Signif_diff_ASVs_ITS[,c(32, 1:31)] 

Signif_diff_ASVs_ITS<-cbind(Signif_diff_ASVs_ITS,Signif_diff_data_ITS)
Signif_diff_ASVs_ITS<-Signif_diff_ASVs_ITS[,1:34] 
write.table(Signif_diff_ASVs_ITS,"Signif_diff_ASVs_ITS.txt", sep = "\t", row.names = FALSE, col.names =TRUE, quote =FALSE) 





