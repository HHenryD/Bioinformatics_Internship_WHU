#args[1]: sample name, such as "OE"
#data: result of dds matrix from DESeq2

#load the packages
library(ggthemes)
library(Cairo)
library(ggplot2)

args <- commandArgs(T)
sample <- args[1]
workdir <- args[2]

#sample <- "OE"
#workdir <- "D:/DESeq2"
setwd(workdir)
result_dir <- paste(workdir,"DEgene", sep="/")
data_file <- paste(workdir, "/", sample, "_allgene.tsv",sep="")
data <- read.table(file=data_file, header=T)

if (!file.exists(result_dir)){dir.create(result_dir, recursive=TRUE)}

data <-data[complete.cases(data),]
data <-as.data.frame(data)
rownames(data) <- data$Row.names
data <- data[,-1]
data$name <- gsub("\\.[0-9]+","",rownames(data))

data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >=1,
                            ifelse((data$log2FoldChange) > 1 ,'Up','Down'),'Not'))
  
##Construct the plot object
# with legend
  
Cairo(file=paste(result_dir,"volcan_PNG_300_lengend_dpi.png",sep="/"), 
  type="png", units="in", bg="white",
  width=5.5, height=5, pointsize=12, dpi=300)

plot <- ggplot(data=data, 
   aes(x=log2FoldChange, y =-log10(padj), 
        colour=threshold,fill=threshold)) +
   scale_color_manual(values=c("blue","grey","red")) +
   geom_point(alpha=0.4, size=1.0) +
   #xlim(c(-30, 30)) +
   #ylim(c(0, 30)) +
   theme_bw(base_size = 12, base_family = "Times") +
   geom_vline(xintercept=c(-0.58,0.58),lty=4,col="grey",lwd=0.6)+
   geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
   theme(legend.position="right",
         panel.grid=element_blank(),
         legend.title = element_blank(),
         legend.text= element_text(face="bold", color="black",family = "Times", size=8),
         plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(face="bold", color="black", size=12),
         axis.text.y = element_text(face="bold",  color="black", size=12),
         axis.title.x = element_text(face="bold", color="black", size=12),
         axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)",y="-log10(padj)",title="Volcano picture of DEG")

print(plot)
ggsave(paste(result_dir,"volcan_PNG_300_lengend_dpi.png",sep="/"), width=5.5, height=5, unit="in",dpi=300)
#pdf(paste(name,"volcan_PNG_300_lengend_dpi.png",sep="_"))
dev.off()


