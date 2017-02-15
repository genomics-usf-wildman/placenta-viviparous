library("data.table")

args <- commandArgs(trailingOnly=TRUE)

## combined_fpkm
load(args[1])

## ensembl_all_ortholog_types
load(args[2])

## cuts -- top 100, top 1%

table_orthologs <- function(this.species,top.100=TRUE) {
    top.results <-
        combined.fpkm[species==this.species &
                      !is.na(gene_id),][order(-FPKM)]
    if (top.100) {
        top.results <-
            top.results[1:100,]
    } else {
        top.results <-
            top.results[1:ceiling(nrow(top.results)*0.01),]
    }
    setkey(top.results,"gene_id")
    species.results <-
        merge(top.results,
              ortholog.types,
              all.x=TRUE,
              all.y=FALSE,
              by.x="gene_id",
              by.y="gene_id")[,table(orthology_type,useNA="ifany")]
    species.results <-
        data.table(species.results)
    species.results[is.na(orthology_type),orthology_type:="novel/unknown"]
    species.results[,species:=this.species]
    if (top.100) {
        species.results[,type:="Top 100"]
    } else {
        species.results[,type:="Top 1%"]
    }
    species.results
}


orthologs.by.species <- list()

for (specie in combined.fpkm[,as.character(unique(species))]) {
    orthologs.by.species[[specie]] <-
        table_orthologs(specie)
    orthologs.by.species[[paste0(specie,"2")]] <-
        table_orthologs(specie,top.100=FALSE)

}

orthologs.by.species <-
    rbindlist(orthologs.by.species)
orthologs.by.species[,ratio:=N/sum(N),by=list(species,type)]

save(orthologs.by.species,
     file=args[length(args)])

