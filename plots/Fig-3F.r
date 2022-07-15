args<-commandArgs(T)
library(ggplot2)
library(ggrepel)
library(ggtext)
library('showtext')
showtext_auto(enable = TRUE)
font_add('Arial','arial.ttf')
font.families()
font.families()


setwd("F:/个性化/HT2020-19711/20220701/20220218_barplot_修改/2022222_another/20220707_frequency/20220712_frequency")
#input <- "Mesenchymal_frequency.xls"
#output <- "Mesenchymal_frequency.pdf"
#order_f <- "Mesenchymal.focal_data_by_genes_res_plot2.xls"

input <- "Metabolic_frequency.xls"
output <- "Metabolic_frequency.pdf"
order_f <- "Metabolic.focal_data_by_genes_res_plot2.xls"

#input <- "Basal_frequency.xls"
#output <- "Basal_frequency.pdf"
#order_f <- "Basal.focal_data_by_genes_res_plot2.xls"

#order_data
data<-read.csv(input, sep='\t')
order_data <- read.table(order_f, sep ='\t',header=T, quote="")
#total_gene_chr <- read.table(total_gene_chr, sep ='\t', header = T)
#total_gene_chr <- total_gene_chr[,c('Gene.Symbol',"Cytoband")]
#high_exp_genes <- read.table(high_exp_genes, sep ='\t', header = T)
#high_exp_genes <- high_exp_genes[,c('Gene.Name')]
#high_exp_genes <- total_gene_chr[total_gene_chr$Gene.Symbol %in% high_exp_genes,]
#data[(data['cytoband'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]
#target_high_exp_genes<- high_exp_genes[high_exp_genes$Cytoband %in% data[(data['cytoband'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]$cytoband,]
#high_exp_genes[high_exp_genes$Cytoband %in% data[(data['cytoband'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]$cytoband,]
#write.table(high_exp_genes,high_exp_output, sep ='\t')
#target_high_exp_genes <- target_high_exp_genes[!duplicated(target_high_exp_genes$Cytoband),]
#if length(target_high_exp_genes)
#  for(i in 1:(length(target_high_exp_genes$Gene.Symbol))){
#    print(paste0(target_high_exp_genes[i,'Cytoband'],"(",target_high_exp_genes[i,'Gene.Symbol'],")"))
#    target_high_exp_genes[i,"Cytoband_gene"] <- paste0(target_high_exp_genes[i,'Cytoband'],"(",target_high_exp_genes[i,'Gene.Symbol'],")")
#  }
#high_exp_genes
data[(data['Cytoband'] != ''),]
data[(data['Cytoband'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]
#high_exp_genes[high_exp_genes$Cytoband %in% data[(data['cytoband'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]$cytoband,]
#colnames(target_high_exp_genes)[2] <- 'cytoband'
#apply(target_high_exp_genes,2,combine)
#target_high_exp_genes
#data$Cytoband_gene <- ""
#target_high_exp_genes
#data <- merge(data,target_high_exp_genes,by.x = 'cytoband', by.y = 'cytoband', all = TRUE)
data[is.na(data)] = ""
#data
#order_data
fa<-unique(order_data$name)
fa2<-unique(order_data$Chr)
data$Chr <- factor(data$Chr,levels=fa2)
data$name <- factor(data$name,levels=fa)
data <- data[order(data$name),]
#View(order_data)
#View(data)
#data[(data['Cytoband_gene'] != '') & ((data['Type'] == 'Gain')|(data['Type'] == 'Amplification')),]

pdf(output,width=10,height=2)
ggplot(data,aes(name,Ratio,fill=Type))+ 
  geom_bar(stat="identity",position='stack',width=0.0001)+theme_classic()+ theme(text = element_text(family = 'Arial'))+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())+
  geom_vline(xintercept=c(-3,3), linetype="dotted")+
  scale_fill_manual(values=c(Deletion="#3653A5",Loss ="#B3BBD6",Gain = "#EBAAAC",Amplification="#E21A21"))+ scale_y_continuous(limits = c(-0.75,0.75), breaks = c(-0.5,0.5)) +
  #facet_grid(.~Chr,scales="free_x",space="free")+  
  theme(strip.text.x = element_text(family ='Arial', size=6,color="black",face="bold.italic"),strip.background = element_rect(color="white",fill="white"),panel.spacing.x=unit(0, "cm")) +
  geom_text_repel(data = data[(data['Cytoband'] != '') & ((data['Type'] == 'Amplification')),],aes(name, as.numeric(Ratio_total),label = Cytoband),family ='Arial', size = 3, segment.size = 0.1,nudge_y = 0.5) + 
  geom_text_repel(data = data[(data['Cytoband'] != '') & ((data['Type'] == 'Loss')|(data['Type'] == 'Deletion')),],aes(name,as.numeric(Ratio_total),label = Cytoband), family ='Arial',size = 3, segment.size = 0.1,nudge_y = -0.5)
dev.off()

