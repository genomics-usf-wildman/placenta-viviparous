args <- commandArgs(trailingOnly=TRUE)
library(seqinr)
library(grid)
cyp19a1.aln <- read.alignment(file=args[1],format="clustal")
source(args[2])
pdf(file=args[length(args)],
    width=10,
    height=10)
dotplot.all(cyp19a1.aln,
            wsize=100,
            wstep=100,
            nmatch=100,
            )
dev.off()
