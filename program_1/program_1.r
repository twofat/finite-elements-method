require(ggplot2)
require(xtable)
setwd("E://Finite Element Methods//program_1//program_1")

output = read.table("norm.txt", head=F)
norm = c(rep(c("max","l2","h1"),4))
output[,2] = log(output[,2]/20) / log(2) 
output = cbind(output, norm)
colnames(output) = c("error", "N", "norm")
xtable(output)

p = ggplot(output,aes(x=N,y=log(error)/log(2),color=norm))
p+geom_point(size=3.5)+geom_line(size=0.5)

par(mfrow=c(2,2))
for(n in c(20,40,80,160)){
outreal = read.table("outreal.txt", head=F)
outestimate = read.table(paste("outestimate",n,".txt",sep=""), head=T)
plot((0:1000)/1000, outreal[,1],xlab="x",ylab="y",cex=0.1,col="pink")
points((0:n)/n,outestimate[,1],col="blue",pch=15,cex=0.3)
title(main=paste("problem_2    N=",n,sep=""))
}






