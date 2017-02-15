library("data.table")


load("ensembl_all_orthologs")

resolve_orthology <- function(orthology_type) {
    orthology_type <- gsub("ortholog_","",
                           orthology_type)
    if (all(orthology_type==orthology_type[1])) {
        return(orthology_type[1])
    }
    if (any(orthology_type=="many2many")) {
        return("many2many")
    }
    if (any(orthology_type=="one2many")) {
        return("one2many")
    }
    return(NA)
}


ortholog.types <-
    orthologs[,list(species=primary[1],
                    orthology_type=resolve_orthology(homolog_orthology_type)),
              by=ensembl_gene_id]
setnames(ortholog.types,
         "ensembl_gene_id",
         "gene_id")
setkey(ortholog.types,
       "gene_id")

save(file="ensembl_all_ortholog_types",
     ortholog.types
     )
