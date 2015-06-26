args <- commandArgs(trailingOnly=TRUE)
library(seqinr)
library(grid)
csh1.aln <- read.alignment(file=args[1],format="clustal")
source(args[2])
pdf(file=args[length(args)],
    width=10,
    height=10)
dotplot.all(csh1.aln,
            start=1,
            end=6910,
            wsize=100,
            wstep=100,
            nmatch=100,
            seq=c(1,2,3,4,5,6,7,8,9,14))
dev.off()
