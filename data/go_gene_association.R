library(data.table)
args <- commandArgs(trailingOnly=TRUE)

go.gene.association <- fread(args[1])
go.gene.association[,V18:=NULL]
go.gene.association[,V19:=NULL]
setnames(go.gene.association,
         c("database",
           "database_id",
           "gene_name",
           "contributes",
           "go_id",
           "go_name",
           "go_namespace",
           "go_ref",
           "go_evidence",
           "other_ids",
           "go_namespace_key",
           "other_names",
           "other_short_names",
           "type",
           "taxon",
           "date",
           "source"))
              
save(file=args[length(args)],
     go.gene.association)
