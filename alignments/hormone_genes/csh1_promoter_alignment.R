args <- commandArgs(trailingOnly=TRUE)
library(seqinr)
library(grid)
csh1.aln <- read.alignment(file=args[1],format="clustal")
source(args[2])
pdf(file=args[length(args)],
    width=10,
    height=10)
dotplot.all(csh1.aln,
            start=72000,
            end=80000,
            wsize=200,
            wstep=200,
            nmatch=200,
            seq=c(4,5,6,7,8,9,10,11,12,13))
dev.off()
