library(openxlsx)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
data<-read.xlsx("data.xlsx",sheet=1,sep.names = " ")
p<-ggplot(data,aes(x = data[,1],y = data[,2])) +
	geom_point(size=4)+
	geom_smooth(method = "lm", fullrange = TRUE) +
	facet_wrap(~pro) +
	theme_bw()+
	theme(text=element_text(family="ArialMT"))+
	stat_cor(data = data,method = "pearson",size = 7,label.y = max(data$value)+1) +
	labs(x = "Protein Abundance of Kinase", y = "Abundance of Substrate Phosphosite") +
	theme(plot.title = element_text(hjust = 0.5, size=20))+
	theme(axis.text = element_text(size = 15,colour = "black"),
	axis.title = element_text(size = 18) ,strip.text = element_text(size = 25))+
	theme(plot.margin = unit(rep(3,4),"lines"))


ggsave("CAMK2B-CAMK2B.pdf", height=10, width=10, plot=p)
ggsave("CAMK2B-CAMK2B.png", type="cairo-png", height=10, width=10, plot=p)
