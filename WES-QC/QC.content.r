#*************************************section 1 usage***************************************

# define usage function
usage=function(){
  
  cat("
      discription:
      plot actgn distribution.
      Email:    wangmin@biomarker.com.cn
      Options: 
          infile 	 	character 	the input file [forced]
          outfile  	character 	the filename for output graph [ default: <infile>.png]
          height 	 	integer 	the height of graph,the unit is mm [optional, default: 80]
          width 	 	integer 	the width of graph, the unit is mm [optional, default: 100]
      example:
          Rscript plot_acgtn.R infile=T01.acgtn
      \n")
  q(status=1);
}


getarg=function(usage=function(){cat("No argument! \n")}){
  #-------------------------------------------------------------------------
  #Description: getopt is primarily intended to be used with “Rscript”. get arguments from Rscript.
  #             return a list data structure containing names of the flags
  #author:     wangmin@biomarker.com.cn
  #---------------------------------------------------------------------------
  arg <- commandArgs(T)
  #
  if(length(arg)==0){
    do.call(usage,list())
    return(NULL)
  }
  arg=read.table(text=arg,sep="=",row.names=1,as.is=T,colClasses ="character")
  arg=data.frame(t(arg),stringsAsFactors=F)
  arg=as.list(arg)
  # ##eval character options to R variable
  for(i in names(arg)){arg[[i]]=type.convert(arg[[i]],as.is=T)}
  arg
}

#******************************************section 2 main function********************************
plot_acgtn=function(infile,outfile=paste(infile,".png",sep=""),height=80,width=120){
  library(RColorBrewer)
  data=read.delim(infile,check=F)
  
  pdf("acgtn.pdf",width=8,height=6)
  par(mar=c(3,3,2,0.5),mgp=c(1.5,0.01,0),tck=-0.005, cex=0.7,cex.lab=0.8,cex.main=0.8,
      xaxs="i",font.lab=2,font.axis=2)
  
  matplot(x=c(1:nrow(data)),y=(data[,-1]*100),type="l",lty=1,col=c("#008B45","#CD1076","#4169E1","#FF7F24","#636363","#00FF00"),axes=F,
          xlim=c(0,nrow(data)),ylim=c(0,100),main="Base Distribution",ylab="Percentage",xlab="",lwd=0.4,font.lab=2,las=2)
  
  
  if(nrow(data)>200){
	read1=which(data[,1]==1)[2]-1
	read2=nrow(data)-read1
	read=nrow(data)
    abline(v=read1,col="darkgray",lwd=0.8)
    pos1=seq(0,read1,50)
    pos2=seq(0,read2,50)
    pos1[1]=1
    pos2=pos2+read1
	pos2[1]=read1+1
	pos3=seq(0,read2,50)
    axis(1,labels=c(pos1,"",pos3[-1]),at=c(pos1,pos2),col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5)
    axis(2,labels=seq(0,100,10),at=seq(0,100,10),cex.axis=0.6,col.axis="gray25",lwd=0.5,las=2)
    mtext("Read1",side=1,line=1.5,at=read1/2,font=2,cex=0.6)
    mtext("Read2",side=1,line=1.5,at=(read1+read2/2),font=2,cex=0.6)
  }else{
    pos=seq(0,nrow(data),25)
    pos[1]=1
    axis(1,labels=pos,at=pos,col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5)
    axis(2,col="gray30",cex.axis=0.6,col.axis="gray25",lwd=0.5,las=2)
    mtext("Cycle Number",side=1,line=1.5,at=max(data[,1]),font=2,cex=0.6)
  }
  
  lgd=sub("\\(%\\)","",colnames(data)[-1])
  legend( "topright",lgd,col=c("#008B45","#CD1076","#4169E1","#FF7F24","#636363","#00FF00"),box.lty=1,bty="n",lwd=1,cex=0.5)

  box(col="gray30",lwd=0.5)
  dev.off()
  
  
}


#************************section 3 call main function*********************************************** 
arg=getarg(usage())

if(!all(c("infile")%in% names(arg))) {usage()}

do.call(plot_acgtn,arg)
