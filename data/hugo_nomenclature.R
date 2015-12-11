library("data.table")

args <- commandArgs(trailingOnly=TRUE)

hugo.nomenclature <- fread(args[1])
setnames(hugo.nomenclature,"Ensembl ID(supplied by Ensembl)","gene_id")
hugo.nomenclature <- hugo.nomenclature[gene_id!="",]
setkey(hugo.nomenclature,"gene_id")

save(hugo.nomenclature,
     file=args[length(args)])
