library(data.table)
library(reshape2)
args <- commandArgs(trailingOnly=TRUE)

fpkm.files <- args[-(1:2)]

output.file <- args[1]

homology.file <- args[2]

reads <- list()
for (fpkm.file in fpkm.files) {
    reads[[fpkm.file]] <- fread(fpkm.file)
    reads[[fpkm.file]][,species:=gsub("^.+Table_\\d+_(.+?)_(?:[Gg]enes|From).+$","\\1",fpkm.file)]
    ## these files have some things misnamed or differently named
    rename.cols <- list("status"="FPKM_status",
                        "HUMAN NAME"="Human Name",
                        "Gene Name"="Human Name",
                        "Human ID"="Human Name",
                        "Human Homolog"="Human Name"
                        )
    for (col.name in names(rename.cols)) {
        if (any(col.name %in% colnames(reads[[fpkm.file]])))  {
            setnames(reads[[fpkm.file]],col.name,rename.cols[[col.name]])
        }
    }
    if (any("other comments" %in% colnames(reads[[fpkm.file]]))) {
        reads[[fpkm.file]][,"other comments":=NULL]
    }
}

### from ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt
homology.mouse.human <- fread(homology.file)
setnames(homology.mouse.human,"HomoloGene ID","homologene.id")
setnames(homology.mouse.human,"Common Organism Name","organism.name")
homology.table <-
    data.table(dcast(homology.mouse.human,
                     homologene.id~organism.name,
                     fun.aggregate=function(x){x[1]},value.var="Symbol"))
setnames(homology.table,"mouse, laboratory","mouse")
setkey(homology.table,"mouse")
mouse.human.homolog <- function(mouse.gene) {
    return(homology.table[mouse==mouse.gene,human])
}

combined.fpkm <- rbindlist(reads,fill=TRUE)
combined.fpkm[species=="Mouse","Human Name":=sapply(gene_short_name,mouse.human.homolog)]
combined.fpkm[,V15:=NULL]
combined.fpkm[,V14:=NULL]



save(file=output.file,combined.fpkm)

