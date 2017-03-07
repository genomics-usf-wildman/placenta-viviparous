library("data.table")

args <- c("housekeeping_genes_superset.txt"
          "housekeeping_genes_superset"
          )

args <- commandArgs(trailingOnly=TRUE)



housekeeping_genes_superset <-
    data.table(read.table(args[1]))

save(housekeeping.genes.superset,file=args[length(args)])
