#该代码获取自：https://github.com/Weidong-Chen-Microbial-Ecology/Stochastic-assembly-of-river-microeukaryotes/blob/master/Neutral%20community%20model.r
#在开始部分的 read.csv() 那里将自己的丰度表读入，然后继续运行后面的代码执行计算即可
#中性群落模型（NCM）的运算细节和 Chen等（2019）的原文献一致，您如果需要变更运算参数或者作图细节等，记得修改 R 代码

library(Hmisc)
library(minpack.lm)
library(stats4)

spp<-read.csv('all_16S.txt',head=T,stringsAsFactors=F,row.names=1,sep = "\t")

#Silk
spp<-spp[,1:15]

#Soil
#spp<-spp[,16:30]


spp<-t(spp)


N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  


    
write.csv(p, file = "Silk_p.csv")   #Silk
#write.csv(p, file = "Soil_p.csv")    #Soil


write.csv(freq, file = "Silk_freq.csv")   #Silk
#write.csv(freq, file = "Soil_freq.csv")   #Soil


write.csv(freq.pred, file = "Silk_freq.pred.csv")  #Silk
#write.csv(freq.pred, file = "Soil_freq.pred.csv")   #Soil


#plot
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('lightgray',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#B3DE69'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#FB8072'
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')

grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 

draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

#Silk
grid.text(y=unit(1.3,'npc')-unit(2.5,'lines'),label='Silk', gp=gpar(fontface=2)) 
#Soil
#grid.text(y=unit(1.3,'npc')-unit(2.5,'lines'),label='Soil', gp=gpar(fontface=2)) 


#save image 4*7inches



rm(list=ls())

#Replot the soil group 


