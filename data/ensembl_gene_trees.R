library(biomaRt)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

species <-
    c("hsapiens",
      "dnovemcinctus",
      "ptroglodytes",
      "btaurus",
      "cfamiliaris",
      "lafricana",
      "mmusculus",
      "mdomestica",
      "ecaballus",
      "sscrofa",
      "oaries"
      )

orthologs <- list()

for (i in 1:length(species)) {
    main.species <- species[i]
    mart <- useMart("ENSEMBL_MART_ENSEMBL",
                    dataset=paste0(main.species,"_gene_ensembl"),
                    host="may2015.archive.ensembl.org")

    orths <-
        data.table(getBM(c("ensembl_gene_id",
                           paste0(secondary.species,
                                  paste0("_homolog_",
                                         c("ensembl_gene",
                                           "orthology_type"))
                                  )),
                         values =TRUE,
                         mart = mart,
                         bmHeader=FALSE))
    setnames(orths,
             paste0(secondary.species,
                                  paste0("_homolog_",
                                         c("ensembl_gene",
                                           "orthology_type"))
                    ),
             c("homolog_ensembl_gene_id","homolog_orthology_type"))
    orths <- orths[homolog_ensembl_gene_id!="",]
    orths[,primary:=main.species]
    orths[,secondary:=secondary.species]
    orthologs[[orth.name]] <- orths
}
orthologs <-
    rbindlist(orthologs)

save(orthologs,
     file=args[length(args)])
