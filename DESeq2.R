#load the packages
library("DESeq2")

args <- commandArgs(T)
fc_dir <- args[1]
treat <- args[2]
ctrl <- args[3]
out_dir <- args[4]
geneid_length <- args[5]

# fc_dir <- "D:/DESeq2"
# treat <- "OE"
# ctrl <- "WT"
# out_dir <- "results"
# geneid_length <- 15

if (!file.exists(out_dir)){dir.create(out_dir, recursive=T)}

treat_data <- c(paste(treat, '1', sep='_'), paste(treat, '2', sep='_'))
ctrl_data <- c(paste(ctrl, '1', sep='_'), paste(ctrl, '2', sep='_'))

data <- c(treat_data, ctrl_data)

#read all featurecounts results
setwd(fc_dir)
allgene <- read.table(paste(data[1], '/', "featureCounts.txt",sep=''), header=F)
data2 <- data[-1]
names(allgene) <- c("gene", data[1])

for (name in data2)
{
  fc_table <- read.table(paste(name, '/', "featureCounts.txt",sep=''), header=F)
  names(fc_table) <- c("gene", name)
  allgene <- merge(allgene, fc_table, by="gene", all=T)
}

#treat the feature count data
us_count <- allgene[,-1]
rownames(us_count) <- allgene[,1]
us_count <- round(us_count, digits=0) 

#prepare the data
us_count <- as.matrix(us_count)
condition <- factor(c(treat, treat, ctrl, ctrl),
                    levels=c(treat, ctrl))
coldata <- data.frame(row.names=colnames(us_count), condition)

#deseq2 analyze
dds <- DESeqDataSetFromMatrix(us_count, coldata, design=~condition) #make dds matrix
dds <- DESeq(dds)             #standardization dds
res <- results(dds)           #get result

#save the result of dds matrix, volcan plot needs it.
setwd(paste("../..", out_dir, sep='/'))
write.table(res, file=paste(treat, "_ddsres.tsv", sep=''),
            quote=F, row.names=T, col.names=T, sep="\t")

#get up/down/all gene from result, set the cut off
res <- res[order(res$padj),]
res <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize=T)),
             by="row.names", sort=F)
deseq_res <- data.frame(res)
up_diff <- subset(deseq_res, (padj < 0.05) & (log2FoldChange > 1))
down_diff <- subset(deseq_res, (padj < 0.05) & (log2FoldChange < -1))
sig_result <- subset(deseq_res, (padj < 0.05) & (abs(log2FoldChange) > 1))
all_result <- subset(deseq_res, baseMean != 0)

#get output files
write.table(up_diff, paste(treat,"_upgene.tsv", sep=''),
            quote=F, row.names=F, sep='\t')
write.table(down_diff, paste(treat,"_downgene.tsv", sep=''),
            quote=F, row.names=F, sep='\t')
write.table(all_result, paste(treat,"_allgene.tsv", sep=''),
            quote=F, row.names=F, sep='\t')

#get gene list
nlen <- as.numeric(geneid_length)
write.table(substr(up_diff$Row.names, 1, nlen), 
            paste(treat, "_upgenelist.txt", sep=''),
            quote=F, row.names=F, col.names=F, sep='\t')
write.table(substr(down_diff$Row.names, 1, nlen), 
            paste(treat, "_downgenelist.txt", sep=''),
            quote=F, row.names=F, col.names=F, sep='\t')
write.table(substr(sig_result$Row.names, 1, nlen), 
            paste(treat, "_siggenelist.txt", sep=''),
            quote=F, row.names=F, col.names=F, sep='\t')

#Volcano Plot
library(ggplot2)
up = read.table("OE_upgene.tsv",header = T)
down = read.table("OE_downgene.tsv",header = T)
all = read.table("OE_allgene.tsv",header = T)
threshold <- as.factor((all$log2FoldChange > 1|all$log2FoldChange < -1)&all$padj<0.05)
ggplot(all,aes(x=log2FoldChange,y=-log10(padj),colour=threshold))+xlab("log2FC")+ylab("-log10padj")+geom_point()
