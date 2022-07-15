args <- commandArgs(T)
file=args[1]
#name=as.character(args[2])
#OutName=paste(file,".png",sep="")
dd <- read.delim(file,sep="\t",header=T)
data <- as.numeric(dd[2:3])
label <- (ceiling(c(data[2],data[1]-data[2])*10000/data[1])/100)
label <- paste(c(data[2],data[1]-data[2]),"(",label,"%)",sep="")
#Main <- paste(name,"QC summary",sep=" ")
#png(OutName,width=700,height=600)
pdf("QCsummary.pdf",width=7,height=6)
par(mar=c(4,2,4,2))
pie(c(data[2],data[1]-data[2]),labels=label,col=c("lightgreen","orchid"),border = NA,main="Summary of Clean Reads" )
legend("bottom",c("High quality reads","Low quality reads"),bg = "white",fill=c("lightgreen","orchid"),bty="n")

dev.off()

