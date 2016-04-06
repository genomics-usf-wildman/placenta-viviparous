library("data.table")
args <- c("../combined_fpkm","all_species_mean_fpkm.txt")
args <- commandArgs(trailingOnly=TRUE)
load(args[1])
write.table(file=args[length(args)],
            combined.fpkm[order(species,-mean_fpkm),],
            sep="\t")

