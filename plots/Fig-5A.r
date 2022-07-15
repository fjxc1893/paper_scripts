library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(ggtext)
library('showtext')
showtext_auto(enable = TRUE)
font_add('Arial','arial.ttf')
font.families()
font.families()

setwd("F:/个性化/HT2020-19711/20211201_HT2020-19711_volcano_bar")

table <- read.table('Basalvsothers.txt', sep="\t", header=T, quote="")
table$q.value <- p.adjust(table$p.value, method = 'BH')
table[,'mark'] <- ''
table[table$q.value < 0.05& table$log2FC > 0, 'mark'] <- 'Up'
table[table$q.value < 0.05& table$FC > 2, 'mark'] <- 'Significatn Up'
table[table$q.value < 0.05& table$log2FC < 0, 'mark'] <- 'Down'
table[table$q.value < 0.05& table$FC < 0.5, 'mark'] <- 'Significatn Down'
write.table(table, 'marked_table.xls', sep = '\t', quote = FALSE, row.names = FALSE)
#go <- table[,c('Accession', 'GO_id', 'GO_term')]
#kegg <- table[,c('Accession','pathway', 'pathway_description')]
#na.omit(kegg)
#write.table(na.omit(go), 'go.backgroud.xls', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(na.omit(kegg), 'kegg.backgroud.xls', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

colnames(table)
table
table <- table[,c('Gene.Name','q.value','log2FC')]
min(table$log2FC)
table['ln_q.value']<- -log10(table$q.value)
table[table['ln_q.value'] >= 12, 'ln_q.value'] <- 12
table
#table[table$log2FC >= 2, 'log2FC'] <- 2
table['group'] <- 'Non-significant'
table[table$q.value < 0.05 & table$log2FC > 0, 'group'] <- 'Up'
table[table$q.value < 0.05 & table$log2FC < 0, 'group'] <- 'Down'
table[table$q.value < 0.05 & table$log2FC > 1, 'group'] <- 'Significant Up'
table[table$q.value < 0.05 & table$log2FC < -1, 'group'] <- 'Significant Down'
table

lfc_up <- length(table[table$group == 'Significant Up','group'])
lfc_down <- length(table[table$group == 'Significant Down','group'])
sig_up <- length(table[table$group == 'Up' ,'group'])
sig_down <- length(table[table$group == 'Down' ,'group'])

table$group <- factor(table$group,levels = c('Significant Up','Significant Down','Up','Down','Non-significant'))
color <- c('#FE0002','#0038FF','#FEBFCA','#86CDF9','#BDBDBD')
p <- ggplot(table, aes(x = log2FC, y = ln_q.value, color = group)) + geom_point() + scale_color_manual(values = color) + geom_vline(xintercept = c(-1,1), linetype = 'longdash', color = 'grey') + geom_hline(yintercept = -log10(0.05), linetype = 'longdash', color = 'grey') + 
  ylim(c(-1,14)) + xlab(expression(paste(log[2]," Fold"," Change",sep =''))) + ylab(expression(paste(-log[10], " q-value", sep =''))) +
  theme(
    panel.border = element_rect(colour = 'black', size = 2,fill = NA),
    panel.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(family = 'Arial', size = 15),
    axis.title = element_text(family = 'Arial', size = 15),
    legend.text = element_text(family = 'Arial'),
    legend.title = element_text(size = 15, family = 'Arial'),
    legend.position = c(0.894,0.131),
    legend.box.background = element_rect(colour = 'black'),
    legend.background = element_blank()
  ) +
  annotate(geom = "richtext",x = c(-2.2,-2.2,-2.2), y = c(4,3.5,3), 
           label = c("<b>Down</b>",paste0("q < 0.05 : ",sig_down),paste0("FC < 0.5 : ",lfc_down)), hjust = 0, fill = NA,label.color = NA, family = 'Arial') + 
  annotate(geom = "richtext",x = c(2,2,2), y = c(4,3.5,3), 
           label = c("<b>Up</b>",paste0("q < 0.05 : ",sig_up),paste0("FC > 2 : ",lfc_up)), hjust = 0, fill = NA,label.color = NA, family = 'Arial')
ggsave("volcano_test.pdf",width = 7.5, height = 7)  

  
