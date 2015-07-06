args <- commandArgs(trailingOnly=TRUE)
library(seqinr)
library(grid)
csh1.aln <- read.alignment(file=args[1],format="fasta")
source(args[2])
pdf(file=args[length(args)],
    width=10,
    height=10)
dotplot.all(csh1.aln,
            wsize=100,
            wstep=1,
            nmatch=100,
            seq=c(4:7,10:11,13:14,3))
dev.off()
