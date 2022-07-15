allargs <- commandArgs(trailingOnly=T)
file <- allargs[1]
#OutName <- allargs[2]
options(scipen = 200)
dd <- read.table(file,sep="\t",header=T)
#png(OutName,width=700,height=600)
pdf("QCsummary.pdf",width=11,height=9)
par(mar=c(4,2,4,2))
# layout(c(1,2),height=c(20,1))
pie(dd[,2],labels=paste(dd[,1],"  (",substring(dd[,3],1,5),"% )",sep=""),col=c("lightgreen","orchid","lightyellow","brown"),border = NA,main="Summary of Clean Reads")
legend("bottom",paste( as.character(dd[,1]),"  (",dd[,2]," , ", substring(dd[,3],1,5),"% )",sep=""),bg = "white",fill=c("lightgreen","orchid","lightyellow","brown"),bty="n",horiz =F)

dev.off()
