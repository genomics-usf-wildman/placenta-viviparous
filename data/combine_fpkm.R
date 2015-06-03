library(data.table)
library(reshape2)
args <- commandArgs(trailingOnly=TRUE)

fpkm.files <- args[-(1:3)]

output.file <- args[1]

homology.file <- args[2]

reads <- list()
for (fpkm.file in fpkm.files) {
    reads[[fpkm.file]] <- fread(fpkm.file)
    reads[[fpkm.file]][,species:=gsub("^.+Table_\\d+_(.+?)_(?:[Gg]enes|From|NEW).+$","\\1",fpkm.file)]
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

### from ftp://ftp.informatics.jax.org/pub/reports/HOM_AllOrganism.rpt
homology.mouse.human <- fread(homology.file)
setnames(homology.mouse.human,"HomoloGene ID","homologene.id")
setnames(homology.mouse.human,"Common Organism Name","organism.name")
homology.table <-
    data.table(dcast(homology.mouse.human,
                     homologene.id~organism.name,
                     fun.aggregate=function(x){x[1]},value.var="Symbol"))
setnames(homology.table,
         c("mouse, laboratory",
           "dog, domestic",
           "rhesus macaque",
           "western clawed frog"
           ),
         c("mouse",
           "dog",
           "macaque",
           "frog"
           ))
human.homolog <- function(gene,species="mouse") {
    return(homology.table[eval(parse(text=paste0(species,"==\"",gene,"\""))),
                          human][1])
}

combined.fpkm <- rbindlist(reads,fill=TRUE)
## CSHL1 is not the same as GH1; call GH1 in Pan an ortholog of CSHL1
setnames(combined.fpkm,
         "Human Name",
         "human_name")

oma.groups <- fread(args[3])
oma.groups[,species:=gsub("\\d+","",oma_entry)]
setkey(oma.groups,"ensembl")
oma.groups.human <- oma.groups[species=="HUMAN"]

oma.homolog <- function(ensembl.id){
    combined.fpkm[tracking_id %in% oma.groups.human[group_num==oma.groups[ensembl.id,group_num],
                                                    ensembl],
                  gene_short_name][1]
}

combined.fpkm[gene_short_name=="GH1 (CSHL1)",gene_short_name:="CSHL1"]
combined.fpkm[gene_short_name=="Q95MK7_PANTR",human_name:="CSHL1"]
combined.fpkm[species=="Homo",human_name:=gene_short_name]
combined.fpkm[human_name=="",human_name:=NA]
combined.fpkm[human_name=="-",human_name:=NA]
combined.fpkm[grepl(" ",human_name),
              human_name:=NA]
### Mouse
combined.fpkm[species=="Mouse" & is.na(human_name),
              human_name := sapply(gene_short_name,function(x){human.homolog(x,species="mouse")})]
### Canis (dog)
combined.fpkm[species=="Canis" & is.na(human_name),
              human_name := sapply(gene_short_name,function(x){human.homolog(x,species="dog")})]
### Cow (cattle)
combined.fpkm[species=="Cow" & is.na(human_name),
              human_name := sapply(gene_short_name,function(x){human.homolog(x,species="cattle")})]
### Pan
combined.fpkm[species=="Pan" & is.na(human_name),
              human_name := sapply(gene_short_name,function(x){human.homolog(x,species="chimpanzee")})]

combined.fpkm[is.na(human_name),
              human_name := sapply(gene_id,oma.homolog)]

combined.fpkm[,V15:=NULL]
combined.fpkm[,V14:=NULL]



### this shouldn't be required, but I've triggered some bug in the
### data.table code
combined.fpkm <- data.table(data.frame(combined.fpkm))

save(file=output.file,combined.fpkm)

