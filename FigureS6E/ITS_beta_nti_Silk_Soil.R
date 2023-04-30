#部分代码来自“微生态”公众号的公开分享
#源代码由刘洪提供，胡天龙进行并行修改
#reference:Yongqi_Shao et al,2022,Science of the Total Environment



rm(list=ls())
library(vegan)
library(Hmisc) 

ASV_all <- read.delim("all_ITS.txt", header = TRUE ,row.names=1)


A=ASV_all 
C=A/rowSums(A)
ASV_all_1<-t(C)

ASV_all_1[ASV_all_1>0]<-1

ASV_all_1<-t(ASV_all_1)
ASV_filtered_ITS<-ASV_all[which(rowSums(ASV_all_1)>15),]



summary(colSums(ASV_filtered_ITS))

ASV_filtered_ITS<-t(ASV_filtered_ITS)
head(rowSums(ASV_filtered_ITS))
ASV_filtered_ITS<-as.data.frame(ASV_filtered_ITS)


ASV_Silk_ITS<- ASV_filtered_ITS[1:15,]
ASV_Soil_ITS<- ASV_filtered_ITS[16:30,]


Fasta_ITS<-read.table("ITS_all.otus.representative.fasta")
re_Fasta_ITS<-Fasta_ITS[seq(1,nrow(Fasta_ITS),2),]
re_Fasta_ITS<-as.data.frame(re_Fasta_ITS)
re_Fasta_ITS$V2<-Fasta_ITS[seq(0,nrow(Fasta_ITS),2),]
colnames(re_Fasta_ITS)<- c("ASVs_names","Sequences")
re_Fasta_ITS$ASVs_names<- gsub(">A*","A",re_Fasta_ITS$ASVs_names) 
row.names(re_Fasta_ITS)<- re_Fasta_ITS$ASVs_names


ASV_filtered_ITS<-t(ASV_filtered_ITS)
ASV_filtered_ITS<-as.data.frame(ASV_filtered_ITS)

ASV_filtered_fasta_ITS<-re_Fasta_ITS[match(row.names(ASV_filtered_ITS), row.names(re_Fasta_ITS)),]

print(row.names(ASV_filtered_fasta_ITS)==row.names(ASV_filtered_ITS))

ASV_filtered_fasta_ITS$ASVs_names<- gsub("ASV",">ASV",ASV_filtered_fasta_ITS$ASVs_names)

write.table(ASV_filtered_fasta_ITS, "ASV_filtered_fasta_ITS.fasta",sep =",", row.names =FALSE, col.names =FALSE, quote =FALSE)

filtered_fasta_ITS<-readLines("ASV_filtered_fasta_ITS.fasta",encoding = "UTF-8")

filtered_fasta_ITS<-gsub(",","\n",filtered_fasta_ITS)

writeLines(filtered_fasta_ITS,"filtered_fasta_ITS.fasta")


#in MacOS command line, perform alignment with muscle
#nohup muscle -in filtered_fasta_ITS.fasta -out aligned_filtered_fasta_ITS.fasta &
# in MacOS command line, construct tree with iqtree
#nohup iqtree -s aligned_filtered_fasta_ITS.fasta -m MFP -B 1000 --bnni -cptime 6000 -T AUTO &


library(treeio)
tree_ITS<-read.tree("aligned_filtered_fasta_ITS.fasta.treefile")


#calculate beta_nti
library(picante)
beta_nti <- function(otu_niche,otu_tree,reps,threads){
  library(picante)
  library(ape)
  
  prune_tree<-prune.sample(otu_niche,otu_tree) 
  
  otu_phydist <- cophenetic(prune_tree)
  
  match.phylo.otu = match.phylo.comm(prune_tree, otu_niche)
 
  beta.mntd.weighted = as.matrix(comdistnt(match.phylo.otu$comm,cophenetic(match.phylo.otu$phy),abundance.weighted=T));
  
  
  beta.reps = reps; # number of randomizations
  rand.weighted.bMNTD.comp = NULL
  dim(rand.weighted.bMNTD.comp)
  library(abind)
  arraybind <- function(...){
    abind(...,along = 3,force.array=TRUE)
  }
  
  library(foreach)
  library(doParallel)
  registerDoParallel(cores = threads)
  rand.weighted.bMNTD.comp <- foreach (rep = 1:beta.reps, .combine = "arraybind") %dopar%{
      
    as.matrix(comdistnt(match.phylo.otu$comm,taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F))
  }
  weighted.bNTI = matrix(c(NA),nrow=nrow(match.phylo.otu$comm),ncol=nrow(match.phylo.otu$comm));
  dim(weighted.bNTI);
  
  for (columns in 1:(nrow(match.phylo.otu$comm)-1)) {
    for (rows in (columns+1):nrow(match.phylo.otu$comm)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals")
    }
  }
  rownames(weighted.bNTI) = rownames(match.phylo.otu$comm)
  colnames(weighted.bNTI) = rownames(match.phylo.otu$comm)
  return(weighted.bNTI)
}



#calculate RC_bray value if site which beta-NTI in (-2,2)
raup_crick_abundance = function(spXsite,
                                plot_names_in_col1 = TRUE,
                                classic_metric = FALSE,
                                split_ties = TRUE,
                                reps = 999,
                                set_all_species_equal = FALSE,
                                as.distance.matrix = TRUE,
                                report_similarity = FALSE,
                                threads = 1) {
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if (plot_names_in_col1) {
    row.names(spXsite) <- spXsite[, 1]
    spXsite <- spXsite[, -1]
  }
  ## count number of sites and total species richness across all plots (gamma)
  library(doParallel)
  library(foreach)
  n_sites <- nrow(spXsite)
  gamma <- ncol(spXsite)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results <-
    matrix(
      data = NA,
      nrow = n_sites,
      ncol = n_sites,
      dimnames = list(row.names(spXsite), row.names(spXsite))
    )
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite / max(spXsite)) -> spXsite.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(spXsite.inc, MARGIN = 2, FUN = sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance <- apply(spXsite, MARGIN = 2, FUN = sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for (null.one in 1:(nrow(spXsite) - 1)) {
    for (null.two in (null.one + 1):nrow(spXsite)) {
      null_bray_curtis <- NULL
      registerDoParallel(cores = threads)
      null_bray_curtis <- foreach (i = 1:reps,.combine = "c") %dopar%{
        ##two empty null communities of size gamma:
        com1 <- rep(0, gamma)
        com2 <- rep(0, gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma,
                    sum(spXsite.inc[null.one, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        com1.samp.sp = sample(which(com1 > 0),
                              (sum(spXsite[null.one, ]) - sum(com1)),
                              replace = TRUE,
                              prob = abundance[which(com1 > 0)])
        com1.samp.sp = cbind(com1.samp.sp, 1)
        # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[, 2], com1.samp.sp[, 1], FUN = sum))
        colnames(com1.sp.counts) = 'counts'
        # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
        # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
        # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp', 'com1.sp.counts')
        ##same for com2:
        com2[sample(1:gamma,
                    sum(spXsite.inc[null.two, ]),
                    replace = FALSE,
                    prob = occur)] <- 1
        com2.samp.sp = sample(which(com2 > 0),
                              (sum(spXsite[null.two, ]) - sum(com2)),
                              replace = TRUE,
                              prob = abundance[which(com2 > 0)])
        com2.samp.sp = cbind(com2.samp.sp, 1)
        # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[, 2], com2.samp.sp[, 1], FUN =
                                                sum))
        colnames(com2.sp.counts) = 'counts'
        # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts))
        # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts
        # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp', 'com2.sp.counts')
        null.spXsite = rbind(com1, com2)
        # null.spXsite;
        ##calculate null bray curtis
        null_bray_curtis = vegdist(null.spXsite, method = 'bray')
        #print(c(i,date()))
      }
      # end reps loop
      ## empirically observed bray curtis
      obs.bray = vegdist(spXsite[c(null.one, null.two), ], method = 'bray')
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis == obs.bray)
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis < obs.bray)
      rc = (num_less_than_in_null) / reps
      # rc;
      if (split_ties) {
        rc = ((
          num_less_than_in_null + (num_exact_matching_in_null) / 2
        ) / reps)
      }
      if (!classic_metric) {
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        rc = (rc - .5) * 2
      }
      results[null.two, null.one] = round(rc, digits = 2)
      ##store the metric in the results matrix
      print(c(null.one, null.two, date()))
    }
    ## end null.two loop
  }
  ## end null.one loop
  if (as.distance.matrix) {
    ## return as distance matrix if so desired
    results <- as.dist(results)
  }
  return(results)
}
## end function






#system.time(beta_nti(otu_niche = ASV_Silk[1:200], otu_tree = tree_ITS,1000,10))

#system.time(beta_nti(otu_niche = ASV_Silk[1:200], otu_tree = tree_ITS,1000,1))

#silk group beta-NTI
paired_beta_nti_silk <- beta_nti(otu_niche = ASV_Silk_ITS, otu_tree = tree_ITS,1000,8)
write.table(paired_beta_nti_silk,"paired_beta_nti_silk.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil group beta-NTI
paired_beta_nti_soil <- beta_nti(otu_niche = ASV_Soil_ITS, otu_tree = tree_ITS,1000,8)
write.table(paired_beta_nti_soil,"paired_beta_nti_soil.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


#calculate RC_bray

#system.time(raup_crick_abundance(otu_niche[,1:200],plot_names_in_col1 = FALSE,threads = 10))

#system.time(raup_crick_abundance(otu_niche[,1:200],plot_names_in_col1 = FALSE,threads = 1))

#system.time(raup_crick_abundance(otu_niche[,1:50],plot_names_in_col1 = FALSE,threads = 1))
#Silk group RC_bray
paired_raup_crick_silk<-raup_crick_abundance(ASV_Silk_ITS,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_silk<-as.matrix(paired_raup_crick_silk)
paired_raup_crick_silk<-as.data.frame(paired_raup_crick_silk)
write.table(paired_raup_crick_silk, "paired_raup_crick_silk.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)
#Soil group RC_bray
paired_raup_crick_soil<-raup_crick_abundance(ASV_Soil_ITS,plot_names_in_col1 = FALSE,threads = 8)
paired_raup_crick_soil<-as.matrix(paired_raup_crick_soil)
paired_raup_crick_soil<-as.data.frame(paired_raup_crick_soil)
write.table(paired_raup_crick_soil, "paired_raup_crick_soil.txt", 
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)





#以下根据Yongqi_Shao老师2022年Science of the Total Environment论文附件扩写
#beta_nti and RCbray
library(reshape2)
library(dplyr)
paired_beta_nti_silk<-read.table(file = "paired_beta_nti_silk.txt", sep = '\t', header = TRUE, row.names = 1)
paired_beta_nti_soil<-read.table(file = "paired_beta_nti_soil.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_silk<- read.table(file = "paired_raup_crick_silk.txt", sep = '\t', header = TRUE, row.names = 1)
paired_raup_crick_soil<- read.table(file = "paired_raup_crick_soil.txt", sep = '\t', header = TRUE, row.names = 1)

dist_melt <- function(dist_a) {
  a <- dist_a 
  a <- as.matrix(a)
  a[upper.tri(a)] <- 10000
  diag(a) <- 10000
  betamat <- melt(a)
  betamat <- betamat[!(betamat$value %in% c(10000)), ]
  return(betamat)
}


beta_nti_silk<- dist_melt(paired_beta_nti_silk)
write.table(beta_nti_silk,"beta_nti_silk_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

beta_nti_soil<- dist_melt(paired_beta_nti_soil)
write.table(beta_nti_soil,"beta_nti_soil_melted.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


raup_crick_silk<- dist_melt(paired_raup_crick_silk)
write.table(raup_crick_silk,"raup_crick_silk.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

raup_crick_soil<- dist_melt(paired_raup_crick_soil)
write.table(raup_crick_soil,"raup_crick_soil.txt",quote = FALSE,
        sep = "\t",row.names = TRUE,col.names = TRUE)









#silk plot



#Silk group Beta_NTI and RCbray results combine

Silk_null_mat <- cbind(beta_nti_silk, raup_crick_silk[, c(3, 2)])[, 1:4]

Silk_null_mat$assembly[Silk_null_mat$value <= -2] <- "Homo Selection"
Silk_null_mat$assembly[Silk_null_mat$value >= 2] <- "Hetro Selection"
Silk_null_mat$assembly[(abs(Silk_null_mat$value) < 2) &
                         Silk_null_mat$value.1 >= 0.95] <-
  "Dispersal Limitation"
Silk_null_mat$assembly[(abs(Silk_null_mat$value) < 2) &
                         Silk_null_mat$value.1 <= -0.95] <-
  "Homogenizing Dispersal"
Silk_null_mat$assembly[(abs(Silk_null_mat$value) < 2) &
                         Silk_null_mat$value.1 > -0.95 &
                         Silk_null_mat$value.1 < 0.95] <- "Drift"
colnames(Silk_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly")

write.table(Silk_null_mat,"Silk_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_Silk_null_mat<- Silk_null_mat %>% group_by(assembly) %>% summarise(n = n())
#Silk_null_mat$assembly






#Silk group Donut plot 
library(ggplot2)
library(RColorBrewer)
pallette_Donut<- brewer.pal(5,'Set1')
# Create data
Silk_Donut_data <- as.data.frame(Summarise_Silk_null_mat)
# Compute percentages
Silk_Donut_data$fraction <- Silk_Donut_data$n / sum(Silk_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
Silk_Donut_data$ymax <- cumsum(Silk_Donut_data$fraction)
# Compute the bottom of each rectangle
Silk_Donut_data$ymin <- c(0, head(Silk_Donut_data$ymax, n=-1))
# Compute label position
Silk_Donut_data$labelPosition <- (Silk_Donut_data$ymax + Silk_Donut_data$ymin) / 2
# Compute a good label
Silk_Donut_data$label <- paste0(Silk_Donut_data$assembly, "\n ", round(Silk_Donut_data$n / sum(Silk_Donut_data$n) *100, 2), "%")


# Make the plot
p1<- ggplot(Silk_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Silk"), size=8) +
  scale_fill_manual(values= c(pallette_Donut)) +
  scale_color_manual(values=c(pallette_Donut)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p1
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI








#Soil group Beta_NTI and RCbray results combine

Soil_null_mat <- cbind(beta_nti_soil, raup_crick_soil[, c(3, 2)])[, 1:4]

Soil_null_mat$assembly[Soil_null_mat$value <= -2] <- "Homo Selection"
Soil_null_mat$assembly[Soil_null_mat$value >= 2] <- "Hetro Selection"
Soil_null_mat$assembly[(abs(Soil_null_mat$value) < 2) &
                         Soil_null_mat$value.1 >= 0.95] <-
  "Dispersal Limitation"
Soil_null_mat$assembly[(abs(Soil_null_mat$value) < 2) &
                         Soil_null_mat$value.1 <= -0.95] <-
  "Homogenizing Dispersal"
Soil_null_mat$assembly[(abs(Soil_null_mat$value) < 2) &
                         Soil_null_mat$value.1 > -0.95 &
                         Soil_null_mat$value.1 < 0.95] <- "Drift"
colnames(Soil_null_mat)=c("Var1","Var2","beta_NTI","RCbray","assembly") 

write.table(Soil_null_mat,"Soil_null_mat.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


Summarise_Soil_null_mat<- Soil_null_mat %>% group_by(assembly) %>% summarise(n = n())
#Soil_null_mat$assembly


Summarise_Soil_null_mat<-as.data.frame(t(Summarise_Soil_null_mat))
Summarise_Soil_null_mat$V5<-c("Hetro Selection", 0)
Summarise_Soil_null_mat<-as.data.frame(t(Summarise_Soil_null_mat))
Summarise_Soil_null_mat<-Summarise_Soil_null_mat[order(Summarise_Soil_null_mat[,1]),] 
Summarise_Soil_null_mat$n<-as.numeric(Summarise_Soil_null_mat$n)




#Soil group Donut plot
library(ggplot2)
library(RColorBrewer)
pallette_Donut<- brewer.pal(5,'Set1')
# Create data
Soil_Donut_data <- as.data.frame(Summarise_Soil_null_mat)
# Compute percentages
Soil_Donut_data$fraction <- Soil_Donut_data$n / sum(Soil_Donut_data$n)
# Compute the cumulative percentages (top of each rectangle)
Soil_Donut_data$ymax <- cumsum(Soil_Donut_data$fraction)
# Compute the bottom of each rectangle
Soil_Donut_data$ymin <- c(0, head(Soil_Donut_data$ymax, n=-1))
# Compute label position
Soil_Donut_data$labelPosition <- (Soil_Donut_data$ymax + Soil_Donut_data$ymin) / 2
# Compute a good label
Soil_Donut_data$label <- paste0(Soil_Donut_data$assembly, "\n ", round(Soil_Donut_data$n / sum(Soil_Donut_data$n) *100, 2), "%")


# Make the plot
p2<- ggplot(Soil_Donut_data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=assembly)) +
  geom_rect() +
  geom_text( x=5, aes(y=labelPosition, label=label, color=assembly), size=6) + # x here controls label position (inner / outer)
  geom_text( x=-1, y=0, label=paste("Soil"), size=8) +
  scale_fill_manual(values= c(pallette_Donut)) +
  scale_color_manual(values=c(pallette_Donut)) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")
p2
#save image 5 inches * 5 inches
#then rearrange the positions of texts with AI






#Merge two pics

library(ggpubr)
p3<-ggarrange(p1,p2, ncol = 2, nrow = 1,common.legend = TRUE,legend="bottom",align = "h")
p3

#save image 6 inches * 10 inches
#then rearrange the positions of texts with AI

