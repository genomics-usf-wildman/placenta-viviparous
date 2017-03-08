library(biomaRt)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="may2015.archive.ensembl.org")


ensembl.orthologs.long <- list()

for (species in c("dnovemcinctus",
                  "ptroglodytes",
                  "btaurus",
                  "cfamiliaris",
                  "lafricana",
                  "mmusculus",
                  "mdomestica")) {
    ensembl.orthologs.long[[species]] <-
        data.table(getBM(c("ensembl_gene_id",
                           paste0(species,"_",
                                  "homolog_ensembl_gene")),
                         values =TRUE, mart = mart, bmHeader=FALSE))
    setnames(ensembl.orthologs.long[[species]],
             paste0(species,"_",
                    "homolog_ensembl_gene"),
             "homolog_ensembl_gene")
    ensembl.orthologs.long[[species]][,species:=species]
    ensembl.orthologs.long[[species]] <-
        ensembl.orthologs.long[[species]][homolog_ensembl_gene!="",]
}

ensembl.orthologs.long <- rbindlist(ensembl.orthologs.long)

save(ensembl.orthologs.long,
     file=args[length(args)])
