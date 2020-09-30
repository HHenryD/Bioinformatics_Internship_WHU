setwd('D:/DESeq2/GO')
up_GO <- read.delim("D:/DESeq2/GO/up_GO.txt")
up_rt = up_GO[up_GO$PValue < 0.05,]
library(tidyr)
up_rt = separate(up_rt, Term, sep = "~",
                 into = c("ID", "Term"))

bp_df = up_rt[up_rt$Category == 'GOTERM_BP_DIRECT',]
bp_df = bp_df[order(bp_df$Count,decreasing = T),]
bp = bp_df[1:10,]

cc_df = up_rt[up_rt$Category == 'GOTERM_CC_DIRECT',]
cc_df = cc_df[order(cc_df$Count,decreasing = T),]
cc = cc_df[1:10,]

mf_df = up_rt[up_rt$Category == 'GOTERM_MF_DIRECT',]
mf_df = mf_df[order(mf_df$Count,decreasing = T),]
mf = mf_df[1:10,]

allGo = rbind(bp,cc,mf)
library(stringr)
table(allGo$Category)
allGo$Category = substr(allGo$Category,8,9)
library(GOplot)
OE_upgene <- read.delim("D:/DESeq2/OE_upgene.tsv")
OE_upgene$Row.names <- sub("(.*)\\..","\\1",OE_upgene$Row.names)
upSig <- OE_upgene[,c(1,3)]
colnames(upSig) <- c('ID','logFC')
data <- allGo[,c(1,2,3,7,6)]
data <- na.omit(data)
colnames(data) = c('category', 'ID', 'term','genes','adj_pval')
circ <- circle_dat(data,upSig)
process<-data$term
chord <- chord_dat(circ, upSig,process)
up_circleplot = GOChord(chord, gene.size = 7, process.label = 10)
ggsave(plot = up_circleplot,'up_circleplot_.png',width = 20,height = 20)
library(ggpubr)
colnames(allGo)
p = ggbarplot(data = allGo,x = "ID",y = 'Count',
              fill = "Category",
              palette = c("cadetblue3","mediumslateblue","mediumorchid3"),
              sort.by.groups = T,xlab = '',ylab = "Target genes")  
ggpar(p,x.text.angle = 90)
ggsave(plot = p,'barplot.png',width = 10,height = 5)
up_circel = GOCircle(circ)
ggsave(plot = up_circle,'up_circle.png',width = 20,height = 20)