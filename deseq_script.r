library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
countdata <- read.table(args[1], header=T, stringsAsFactors=F)
genenames <- countdata$Geneid
countdata <- countdata[, 7:ncol(countdata)] 
colnames(countdata) <- c(paste("G1_rep_",1:3,sep=""), paste("G2_rep_",1:3,sep="")) #rename columns
countdata <- as.matrix(countdata)
rownames(countdata) <- genenames

#Get a boolean vector to remove all ERCC entries
sel <- sapply(rownames(countdata), function(x){ if(substr(x, 1,5)=="ERCC-"){return(FALSE)}else{return(TRUE)} })
countdata <- countdata[sel, ]

coldata <- data.frame("condition"=as.factor(c(rep("cancer", 3), rep("ref", 3))), row.names=colnames(countdata))

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
write.table(res, file=args[2], col.names=T, row.name=F, quote=F)
