args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(reshape)
library(ggplot2)

load(args[1])

library(biomaRt)

mart <- useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

results <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol",
                     "go_id","name_1006","go_linkage_type","namespace_1003"),
                 filters="ensembl_gene_id",
                 values=combined.fpkm[species=="Homo",tracking_id],
                 mart=mart)
go.functional.groups <- data.table(results)

save(file=args[length(args)],
     go.functional.groups)
