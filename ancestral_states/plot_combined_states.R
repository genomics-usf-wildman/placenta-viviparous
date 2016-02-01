library("data.table")
library("phytools")
library("ape")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])

pdf(file=args[length(args)])

for(gene in names(combined.anc.states)) {
    try({contMap(combined.anc.states[[gene]]$tree,
                 combined.anc.states[[gene]]$expression.vector)
                 title(gene)
    },silent=TRUE)
}

     
dev.off()
