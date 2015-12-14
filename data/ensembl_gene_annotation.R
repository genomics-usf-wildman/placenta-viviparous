library("biomaRt")
library("data.table")

args <- commandArgs(trailingOnly=TRUE)



gene.annotation <- list()
for (species in c("hsapiens","btaurus","cfamiliaris","dnovemcinctus","ecaballus",
                  "lafricana","mdomestica","mmusculus","oaries","ptroglodytes",
                  "sscrofa")) {
    mart <- useMart("ensembl",dataset=paste0(species,"_gene_ensembl"))

    gene.annotation[[species]] <-
        data.table(getBM(c("ensembl_gene_id",
                           "external_gene_name",
                           "gene_biotype",
                           "chromosome_name",
                           "start_position"
                           ),
                         value=TRUE,mart=mart,bmHeader=FALSE))
    
    gene.annotation[[species]][,species:=species]
}
gene.annotation <- rbindlist(gene.annotation)
save(gene.annotation,
     file=args[length(args)])
