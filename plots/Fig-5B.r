library(ggplot2)
library(ggtext)
library('showtext')
showtext_auto(enable = TRUE)
font_add('Arial','arial.ttf')
font.families()
font.families()

setwd("C:/Users/oebio/Desktop/个性化/HT2020-19711/20211201_HT2020-19711_volcano_bar")

table <- read.table('targetPath.xls', sep="\t", header=T, quote="")
table['Enrichment_score'] <- (table$S_Gene_Number/table$TS_Gene_Number)/(table$B_Gene_Number/table$TB_Gene_Number)
table
ggplot(table, aes(x = Enrichment_score, y = reorder(Pathway_Name,Enrichment_score), fill = p_Value)) + geom_bar(stat = 'identity', position = 'dodge')  + scale_fill_gradientn(colours = c('#4575B5', '#E0F3F7',  '#FFF6B1', 'red')) + ylab("") + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(family = 'Arial'),
        legend.text = element_text(family = 'Arial'),
        legend.title = element_text(family = 'Arial'),
        axis.title = element_text(family = 'Arial'))
pdf(file = "bar_plot_2.pdf",  height=2.92,width=8.59, family = 'sans')
ggplot(table, aes(x = Enrichment_score, y = reorder(Pathway_Name,Enrichment_score), fill = p_Value)) + geom_bar(stat = 'identity', position = 'dodge', width = 0.85)  + scale_fill_gradientn(colours = c('#FFE4E1','#FF3030','#A52A2A')) + xlab("Enrichment score") +ylab("") + guides(fill = guide_colorbar(title = "pvalue")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(family = 'Arial', size = 12, color = 'black'),
        legend.text = element_text(family = 'Arial', color = 'black'),
        legend.title = element_text(family = 'Arial', color = 'black'),
        axis.title = element_text(family = 'Arial', size = 12, color = 'black'),
        axis.title.x = element_text(family = 'Arial', size = 15, color = 'black'))

dev.off()

pdf(file = "bar_plot_3.pdf",  height=2.92,width=8.59, family = 'sans')
ggplot(table, aes(x = Enrichment_score, y = reorder(Pathway_Name,Enrichment_score), fill = p_Value)) + geom_bar(stat = 'identity', position = 'dodge', width = 0.85)  + scale_fill_gradient(low ='#FFE4E1', high = '#A52A2A') + xlab("Enrichment score") +ylab("") + guides(fill = guide_colorbar(title = "pvalue")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text = element_text(family = 'Arial', size = 12, color = 'black'),
        legend.text = element_text(family = 'Arial', color = 'black'),
        legend.title = element_text(family = 'Arial', color = 'black'),
        axis.title = element_text(family = 'Arial', size = 12, color = 'black'),
        axis.title.x = element_text(family = 'Arial', size = 15, color = 'black'))

dev.off()
