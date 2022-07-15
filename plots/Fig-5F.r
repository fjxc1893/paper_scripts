library("ggplot2")
setwd("F:/个性化/HT2020-19711/20211109_complexPlot")

windowsFonts(Arial=windowsFont("Arial"))

#diff Table
group <- read.table(normalizePath("diff_group.txt"), header = T, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "", colClasses = "character")
bassal_s <- unlist(strsplit(as.character(group$control), ','))
other_s <- unlist(strsplit(as.character(group$case), ','))
#control_name <- c(as.character(group$control_name))
#case_name <- c(as.character(group$case_name))
samples <- c(bassal_s, other_s)
bassal_s <- gsub("-",".",bassal_s)
other_s <- gsub("-",".",other_s)

#Q1 target Genes
B_markers <- c('CHGA','CHGB','SCG2','SCG3','SCGN','NCAM1','SYP','SNAP25','UCHL1')
MTC_markers <- c('CALCB','CEACAM5','CEACAM6')
#proocessing Pro
pre_pro <- read.table("post_pro.txt", sep = '\t', header = TRUE,quote = "")
 pre_pro['q.value'] <- p.adjust(pre_pro$p.value, method = "BH")
pre_pro
pro_t <- pre_pro[which(pre_pro$Gene.Name %in% B_markers|pre_pro$Gene.Name %in% MTC_markers),]
write.table(pro_t,'Pro.xls', sep ='\t', quote = FALSE,row.names = FALSE)
B_pro <- pre_pro[which(pre_pro$Gene.Name %in% B_markers),][,c('Gene.Name','log2FC','p.value','q.value')]
#B_pro$q.value <- p.adjust(B_pro$p.value, method = "BH")
rownames(B_pro) <- B_pro[,'Gene.Name']
B_pro <- B_pro[B_markers,]

pre_pro[which(pre_pro$Gene.Name == 'SCGN'),]
M_pro <- pre_pro[which(pre_pro$Gene.Name %in% MTC_markers),][,c('Gene.Name','log2FC','p.value','q.value')]
#M_pro$q.value <- p.adjust(M_pro$p.value, method = "BH")
rownames(M_pro) <- M_pro[,'Gene.Name']
M_pro <- M_pro[MTC_markers,]
colnames(B_pro) <- c('gene_id','log2FoldChange', 'p.value','q.value')
colnames(M_pro) <- c('gene_id','log2FoldChange', 'p.value','q.value')
B_pro['group1'] <- 'Broad spectrum neuroendocrine markers'
M_pro['group1'] <- 'MTC-specific markers'

pro <- rbind(B_pro,M_pro)
pro['group2'] <- 'protein'

# mRNA
pre_mRNA <- read.table("Basal-vs-Other.xls", sep = '\t', header = TRUE)
colnames(pre_mRNA)
mRNA_t <- pre_mRNA[which(pre_mRNA$gene_id %in% B_markers|pre_mRNA$gene_id %in% MTC_markers),]
write.table(mRNA_t,'mRNA.xls', sep ='\t', quote = FALSE,row.names = FALSE)
mRNA <- pre_mRNA[,c('gene_id','log2FoldChange','p.value', 'q.value')]
mRNA1 <- mRNA[which(mRNA$gene_id %in% B_markers),]
mRNA1[which(mRNA1$gene_id %in% B_markers),'group1'] <- 'Broad spectrum neuroendocrine markers'
rownames(mRNA1) <- mRNA1[,'gene_id']
mRNA1 <- mRNA1[B_markers,]

mRNA2 <- mRNA[which(mRNA$gene_id %in% MTC_markers),]
mRNA2[which(mRNA2$gene_id %in% MTC_markers),'group1'] <- 'MTC-specific markers'
rownames(mRNA2) <- mRNA2[,'gene_id']
mRNA2 <- mRNA2[MTC_markers,]
mRNA <- rbind(mRNA1,mRNA2)
mRNA['group2'] <- 'mRNA'
mRNA['CHGB','log2FoldChange'] <- 0
mRNA['CEACAM6','log2FoldChange'] <- 0

#plot log2FC
data <- rbind(mRNA,pro)
tile_mean <- mean(data$log2FoldChange)
lowest_v <- min(data$log2FoldChange)
largest_v <- max(data$log2FoldChange)
mRNA1 <- na.omit(mRNA[which(data['group1'] == 'Broad spectrum neuroendocrine markers'),])
mRNA2 <- na.omit(mRNA[which(data['group1'] == 'MTC-specific markers'),])
pre_mRNA
ggplot(mRNA, aes(x = gene_id,y = log2FoldChange, fill = log2FoldChange)) + geom_bar(stat = 'identity', position = 'stack')
library('reshape2')
test <- pre_mRNA[,c('gene_id',"log2FoldChange",'p.value','q.value')]
test <- melt(test)
ggplot(test,aes(x = variable, y = gene_id, fill = value)) + geom_tile()
colnames(pre_mRNA)

library(forcats) 
rownames(mRNA1) <- 1:nrow(mRNA1)
mRNA_p1 <- ggplot(mRNA1,aes(1,y = factor(gene_id,levels = rev(gene_id)), fill = log2FoldChange)) + geom_tile(color = 'white') +  theme_bw() + scale_fill_gradient(low = '#FEB69E', high = 'red', limits = c(lowest_v,largest_v))+ guides(fill = "none") + 
  theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+   coord_fixed() 

        #axis.text.x = element_blank(),
        #axis.text.y = element_blank())
mRNA_p1
mRNA_p2 <- ggplot(mRNA2,aes(1,y = factor(gene_id,levels = rev(gene_id)), fill = log2FoldChange)) + geom_tile(color = 'white') +   coord_fixed() +
  #facet_wrap(group1, strip.position="left")+
  theme_bw() + scale_fill_gradient(low = '#FEB69E', high = 'red', limits = c(lowest_v,largest_v))+ guides(fill = "none") + 
  theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

#axis.text.x = element_blank(),
#axis.text.y = element_blank())


pro1 <- na.omit(pro[which(pro['group1'] == 'Broad spectrum neuroendocrine markers'),])
pro2 <- na.omit(pro[which(pro['group1'] == 'MTC-specific markers'),])


pro_p1 <- ggplot(pro1,aes(1,y = factor(gene_id,levels = rev(gene_id)), fill = log2FoldChange)) + geom_tile(color = 'white') +   coord_fixed() +  scale_fill_gradient(low = '#FEB69E', high = 'red', limits = c(lowest_v,largest_v))+ guides(fill = "none") + theme_bw() + 
  theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',family = 'Arial', size =13)) + scale_y_discrete(position = "right") 


pro_p1
pro_p2 <- ggplot(pro2,aes(1,y = factor(gene_id,levels = rev(gene_id)), fill = log2FoldChange)) + geom_tile(color = 'white') +   coord_fixed() +  theme_bw() + scale_fill_gradient(low = '#FEB69E', high = 'red', limits = c(lowest_v,largest_v))+ guides(fill = "none") + 
  theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',family = 'Arial', size =13)
        ) + scale_y_discrete(position = "right")

pro_p2

# plot bar plot
total <- rbind(mRNA,pro)
total['log10qvalue'] <- -log(as.numeric(total$q.value), 10)
B_total <- total[which(total['group1'] == 'Broad spectrum neuroendocrine markers'),]
M_total <- total[which(total['group1'] == 'MTC-specific markers'),]
B_total
total_p1 <-  ggplot(data =B_total, aes(x = as.numeric(log10qvalue), y = factor(gene_id,levels = rev(B_markers)), fill = group2)) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') +theme_bw() + guides(fill = 'none') + xlab(expression(paste(-log[10],qvalue, sep =''))) +  scale_fill_manual(values=c('#A8CF98','#FFDE7F')) + theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(size = 1.5,color = 'black'),
        axis.line.x = element_line(size = 1.5,color = 'black'),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 15, color = 'black', vjust = 25,family = 'Arial'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15, color = 'black',family = 'Arial')
        )  + scale_x_continuous(expand=c(0,0), position = 'top',limits = c(0,10)) + geom_vline(aes(xintercept = -log(0.05,10)),size = 1,linetype = 'dashed')
total_p1
#total_p1
#M_total
#B_total
B_total
M_total
MTC_markers
factor(M_total$gene_id,levels = MTC_markers)
MTC_markers
total_p2 <-  ggplot(data =M_total, aes(x = as.numeric(log10qvalue), y = factor(gene_id,levels = c("CEACAM6","CEACAM5","CALCB")), fill = group2 )) + geom_bar(stat = 'identity', position = 'dodge', color = 'black') +theme_bw() + guides(fill = 'none') + 
  theme(panel.border = element_blank(),
        panel.background =element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(size = 1.5,color = 'black'),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) + 
  scale_fill_manual(values=c('#F4AF92','#5D9AD1'))  +scale_x_continuous(expand=c(0,0),limits = c(0,10)) + geom_vline(aes(xintercept = -log(0.05,10)),size = 1,linetype = 'dashed') 
total_p2
#total$log2FoldChange
#lowest_v
total_p3 <-total_p2 + geom_point(data = total, aes(color =log2FoldChange , shape = NA)) + scale_colour_gradient(low = '#FEB69E', high = 'red', limits = c(lowest_v,largest_v),breaks = c(lowest_v,largest_v), label = c('0','1.8'),guide = 'colorbar',position = 'bottom') + theme(legend.position = 'bottom',legend.title = element_blank(),legend.text = element_text(family = 'Arial',size = 12))
#  guides(fill = guide_colorbar(ticks = FALSE,direction = 'horizontal'))
total_p3

#total
library('patchwork')
#library('ggpubr')

#ggarrange(mRNA_p1,mRNA_p2,ncol = 1)

library('showtext')
showtext_auto(enable = TRUE)
font_add('Arial','arial.ttf')
font.families()
font.families()
pdf(file = "plot1.pdf",  height=8,width=12, family = 'sans')
((mRNA_p1/mRNA_p2)|(pro_p1/pro_p2))|(total_p1/total_p2 + plot_layout(heights = c(9,3)))
dev.off()

tiff(filename = "plot1.tiff", res = 400,units = 'px',  height=16,width=18)
((mRNA_p1/mRNA_p2)|(pro_p1/pro_p2))|(total_p1/total_p2 + plot_layout(heights = c(9,4)))
dev.off()

#total_p/(geom_point(data = total, aes(size="log10pvalue", shape = NA), colour = "grey50"))
write.table(total,'total.xls', sep ='\t', quote = FALSE,row.names = FALSE)
total




  #+ plot_layout(widths = c(2, 1))



#+ scale_y_discrete(position = "right") + scale_x_discrete(position = "top")
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #      axis.ticks=element_blank(), axis.text=element_text(size=5), plot.title=element_text(hjust=0),
  #      strip.text=element_text(hjust=0), panel.margin.x=unit(0.5, "cm"),panel.margin.y=unit(0.5, "cm"),
  #      )

