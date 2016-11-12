library("data.table")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
write.table(placenta.classification.p[fdr<0.05 & ((type=="aov" & factor=="shape")| type=="polr"),unique(gene)],
            quote=FALSE,
            sep="",
            col.names=FALSE,
            row.names=FALSE,
            file=args[length(args)])
