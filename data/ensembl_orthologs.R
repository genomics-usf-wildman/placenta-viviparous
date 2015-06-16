library(biomaRt)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

human <- useMart("ensembl",dataset="hsapiens_gene_ensembl")


orthologs <-
    data.table(getBM(c("ensembl_gene_id",
                       paste(sep="_",
                             c("dnovemcinctus",
                               "ptroglodytes",
                               "btaurus",
                               "cfamiliaris"
                               ),
                             "homolog_ensembl_gene")),
                     values =TRUE, mart = human, bmHeader=FALSE))

orthologs.2 <-
    data.table(getBM(c("ensembl_gene_id",
                       paste(sep="_",
                             c("lafricana",
                               "mmusculus",
                               "mdomestica"),
                             "homolog_ensembl_gene")),
                     values =TRUE, mart = human, bmHeader=FALSE))

orthologs <- orthologs[!duplicated(ensembl_gene_id),]
orthologs.2 <- orthologs.2[!duplicated(ensembl_gene_id),]
setkey(orthologs,"ensembl_gene_id")
setkey(orthologs.2,"ensembl_gene_id")
orthologs <-
    orthologs[orthologs.2]

save(orthologs,
     file=args[length(args)])
