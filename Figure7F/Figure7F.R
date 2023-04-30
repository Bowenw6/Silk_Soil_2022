#来自微生态公众号
#源代码由刘洪提供，胡天龙进行并行修改。
#矩阵数据合并等部分代码根据Yongqi_Shao老师2022年Science of the Total Environment论文附件扩写



rm(list=ls())
library(vegan)
library(Hmisc) 

ASV_all <- read.delim("all_16S.txt", header = TRUE ,row.names=1)

A=ASV_all 
C=A/rowSums(A)
ASV_all_1<-t(C)

ASV_all_1[ASV_all_1>0]<-1

ASV_all_1<-t(ASV_all_1)
ASV_filtered_16S<-ASV_all[which(rowSums(ASV_all_1)>15),]



summary(colSums(ASV_filtered_16S))

ASV_filtered_16S<-t(ASV_filtered_16S)
head(rowSums(ASV_filtered_16S))
ASV_filtered_16S<-as.data.frame(ASV_filtered_16S)

ASV_Silk_16S<- ASV_filtered_16S[1:15,]
ASV_Soil_16S<- ASV_filtered_16S[16:30,]


ASV_Silk_16S_2w<-ASV_Silk_16S[1:3,]
ASV_Silk_16S_4w<-ASV_Silk_16S[4:6,]
ASV_Silk_16S_6w<-ASV_Silk_16S[7:9,]
ASV_Silk_16S_8w<-ASV_Silk_16S[10:12,]

ASV_Soil_16S_2w<-ASV_Soil_16S[1:3,]
ASV_Soil_16S_4w<-ASV_Soil_16S[4:6,]
ASV_Soil_16S_6w<-ASV_Soil_16S[7:9,]
ASV_Soil_16S_8w<-ASV_Soil_16S[10:12,]


library(treeio)
tree_16S<-read.tree("aligned_filtered_fasta_16S.fasta.treefile")


#compute beta_nti
library(picante)
beta_nti <- function(otu_niche,otu_tree,reps,threads){
  library(picante)
  library(ape)
  
  prune_tree<-prune.sample(otu_niche,otu_tree) 
  
  otu_phydist <- cophenetic(prune_tree)
  
  match.phylo.otu = match.phylo.comm(prune_tree, otu_niche)
  
  beta.mntd.weighted = as.matrix(comdistnt(match.phylo.otu$comm,cophenetic(match.phylo.otu$phy),abundance.weighted=T));
  
  #identical(rownames(match.phylo.otu$comm),colnames(beta.mntd.weighted));
  # just a check, should be TRUE
  #identical(rownames(match.phylo.otu$comm),rownames(beta.mntd.weighted)); # just a check, should be TRUE
  # calculate randomized betaMNTD
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






#system.time(beta_nti(otu_niche = ASV_Silk[1:200], otu_tree = tree_16S,1000,10))

#system.time(beta_nti(otu_niche = ASV_Silk[1:200], otu_tree = tree_16S,1000,1))


#silk_2w beta-NTI
paired_beta_nti_silk_2w <- beta_nti(otu_niche = ASV_Silk_16S_2w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_silk_2w,"paired_beta_nti_silk_2w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_4w beta-NTI
paired_beta_nti_silk_4w <- beta_nti(otu_niche = ASV_Silk_16S_4w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_silk_4w,"paired_beta_nti_silk_4w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_6w beta-NTI
paired_beta_nti_silk_6w <- beta_nti(otu_niche = ASV_Silk_16S_6w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_silk_6w,"paired_beta_nti_silk_6w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_8w beta-NTI
paired_beta_nti_silk_8w <- beta_nti(otu_niche = ASV_Silk_16S_8w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_silk_8w,"paired_beta_nti_silk_8w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)



#soil_2w beta-NTI
paired_beta_nti_soil_2w <- beta_nti(otu_niche = ASV_Soil_16S_2w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_soil_2w,"paired_beta_nti_soil_2w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil_4w beta-NTI
paired_beta_nti_soil_4w <- beta_nti(otu_niche = ASV_Soil_16S_4w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_soil_4w,"paired_beta_nti_soil_4w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


#soil_6w beta-NTI
paired_beta_nti_soil_6w <- beta_nti(otu_niche = ASV_Soil_16S_6w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_soil_6w,"paired_beta_nti_soil_6w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)


#soil_8w beta-NTI
paired_beta_nti_soil_8w <- beta_nti(otu_niche = ASV_Soil_16S_8w, otu_tree = tree_16S,1000,8)

write.table(paired_beta_nti_soil_8w,"paired_beta_nti_soil_8w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)






#以下根据Yongqi_Shao老师2022年Science of the Total Environment论文附件扩写
#计算Silk组和Soil组的beta_nti和RCbray
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

#silk_2w_beta_nti
beta_nti_silk_2w<- dist_melt(paired_beta_nti_silk_2w)
write.table(beta_nti_silk_2w,"beta_nti_silk_melted_2w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_4w_beta_nti
beta_nti_silk_4w<- dist_melt(paired_beta_nti_silk_4w)
write.table(beta_nti_silk_4w,"beta_nti_silk_melted_4w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_6w_beta_nti
beta_nti_silk_6w<- dist_melt(paired_beta_nti_silk_6w)
write.table(beta_nti_silk_6w,"beta_nti_silk_melted_6w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#silk_8w_beta_nti
beta_nti_silk_8w<- dist_melt(paired_beta_nti_silk_8w)
write.table(beta_nti_silk_8w,"beta_nti_silk_melted_8w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil_2w_beta_nti
beta_nti_soil_2w<- dist_melt(paired_beta_nti_soil_2w)
write.table(beta_nti_soil_2w,"beta_nti_soil_melted_2w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil_4w_beta_nti
beta_nti_soil_4w<- dist_melt(paired_beta_nti_soil_4w)
write.table(beta_nti_soil_4w,"beta_nti_soil_melted_4w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil_6w_beta_nti
beta_nti_soil_6w<- dist_melt(paired_beta_nti_soil_6w)
write.table(beta_nti_soil_6w,"beta_nti_soil_melted_6w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)

#soil_8w_beta_nti
beta_nti_soil_8w<- dist_melt(paired_beta_nti_soil_8w)
write.table(beta_nti_soil_8w,"beta_nti_soil_melted_8w.txt",quote = FALSE,
            sep = "\t",row.names = TRUE,col.names = TRUE)




#Merge results

beta_nti_silk_all<- rbind(beta_nti_silk_2w, beta_nti_silk_4w, beta_nti_silk_6w, beta_nti_silk_8w)
beta_nti_silk_all$Time<-c(rep("2", 3), rep("4", 3), rep("6", 3), rep("8", 3))
beta_nti_silk_all<-beta_nti_silk_all[,3:4]

colnames(beta_nti_silk_all)<-c("Silk","Weeks")


beta_nti_soil_all<- rbind(beta_nti_soil_2w, beta_nti_soil_4w, beta_nti_soil_6w, beta_nti_soil_8w)
beta_nti_soil_all<-beta_nti_soil_all[,1:3]


beta_nti_lm_plot<- cbind(beta_nti_silk_all, beta_nti_soil_all)

beta_nti_lm_plot<-beta_nti_lm_plot[,c(TRUE,TRUE,FALSE,FALSE,TRUE)]

colnames(beta_nti_lm_plot)<-c("Silk", "Weeks","Soil")

beta_nti_lm_plot<-beta_nti_lm_plot[,c("Silk", "Soil","Weeks")]
rownames(beta_nti_lm_plot)<-c(rep(1:12))

beta_nti_lm_plot<-reshape2::melt(beta_nti_lm_plot, id= 'Weeks')
beta_nti_lm_plot$Weeks<-as.numeric(beta_nti_lm_plot$Weeks)

colnames(beta_nti_lm_plot)<-c("Weeks", "Groups", "BetaNTI")


#lm plot
library(ggplot2)
library(ggpubr)

p1<-ggplot(beta_nti_lm_plot %>% na.omit()) +
  geom_point(aes(Weeks, BetaNTI, colour = Groups, shape = c("2")), 
             size = 2) +
  geom_smooth(aes(Weeks, BetaNTI, colour = Groups), method = lm) +
  stat_cor(aes(Weeks, BetaNTI, colour = Groups), method = "pearson")+
  
  scale_colour_manual(values = c("Silk" = "#55AFFF",
                                 "Soil"="#F4B800"
                                 ))+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank())+ 
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.title.y=element_text(colour='black', size=18,face = "bold",vjust = 1.5),
        axis.text.y=element_text(colour='black',size=15),
        axis.text.x=element_text(colour = "black",size = 15,
                                 hjust = 1,vjust = 0.5))

p1


#save image 4*3 inches







