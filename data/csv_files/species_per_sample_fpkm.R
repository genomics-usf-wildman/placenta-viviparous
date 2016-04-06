library("data.table")
args <- c("../combined_fpkm_per_sample","all_species_per_sample_mean_fpkm.txt")
args <- commandArgs(trailingOnly=TRUE)
load(args[1])
write.table(file=args[length(args)],
            combined.fpkm[order(species,-mean_fpkm,file),],
            sep="\t")

