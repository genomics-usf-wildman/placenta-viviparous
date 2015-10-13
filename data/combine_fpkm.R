library(data.table)
library(reshape2)
library(parallel)

args <- commandArgs(trailingOnly=TRUE)

if (is.null(getOption("mc.cores"))) {
    options(mc.cores=12)
}

output.file <- args[1]

homology.file <- args[2]

load(args[4])


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

oma.groups <- fread(args[3])
oma.groups[,species:=gsub("\\d+","",oma_entry)]
setkey(oma.groups,"ensembl")
oma.groups.human <- oma.groups[species=="HUMAN"]

oma.homolog <- function(ensembl.id){
    gene.fpkms[tracking_id %in% oma.groups.human[group_num==oma.groups[ensembl.id,group_num],
                                                    ensembl],
                  gene_short_name][1]
}

### ensembl_orthologs
load(args[5])

orthologs.long <- melt(orthologs,id.vars="ensembl_gene_id")
orthologs.long[,variable:=NULL]
setkey(orthologs.long,"value")

ensembl.homolog <- function(ensembl.id) {
    hum.id <- orthologs.long[ensembl.id,ensembl_gene_id][1]
    if (is.null(hum.id) || is.na(hum.id))
        return(NA)
    return(gene.fpkms[tracking_id==hum.id,gene_short_name][1])
}

gene.fpkms[species=="homo sapiens",human_name:=gene_short_name]
gene.fpkms[human_name=="",human_name:=NA]
gene.fpkms[human_name=="-",human_name:=NA]
gene.fpkms[grepl(" ",human_name),
              human_name:=NA]
### ### Mouse
### gene.fpkms[species=="mus musculus" & is.na(human_name),
###            human_name := mcmapply(function(x){human.homolog(x,species="mouse")},gene_short_name)]
### ### Canis (dog)
### gene.fpkms[species=="canis familiaris" & is.na(human_name),
###               human_name := mcmapply(function(x){human.homolog(x,species="dog")},gene_short_name)]
### ### Cow (cattle)
### gene.fpkms[species=="bos taurus" & is.na(human_name),
###               human_name := mcmapply(function(x){human.homolog(x,species="cattle")},gene_short_name)]
### ### pan paniscus
### gene.fpkms[species=="pan paniscus" & is.na(human_name),
###               human_name := mcmapply(function(x){human.homolog(x,species="chimpanzee")},gene_short_name)]

gene.fpkms[is.na(human_name),
           human_name := mcmapply(ensembl.homolog,gene_id)]

gene.fpkms[is.na(human_name),
              human_name := mcmapply(oma.homolog,gene_id)]

### this shouldn't be required, but I've triggered some bug in the
### data.table code
combined.fpkm <- data.table(data.frame(gene.fpkms))

setkey(oma.groups,"ensembl")
setkey(combined.fpkm,"tracking_id")
combined.fpkm <- oma.groups[,list(group_num,fingerprint,ensembl)][combined.fpkm]
combined.fpkm[,ensembl:=NULL]
setnames(combined.fpkm,"group_num","oma_group_num")

combined.fpkm[,oma_group_name:=as.character(NA)]
combined.fpkm[!is.na(oma_group_num),
              oma_group_name:=paste(collapse=",",
                  grep("^-$",unique(gene_short_name),value=TRUE,invert=TRUE)),
              by=as.numeric(oma_group_num)]

combined.fpkm[,c("class_code","nearest_ref_id",
                 "i.gene_short_name",
                 "tss_id","locus","length","coverage"):=
                     list(NULL,NULL,NULL,NULL,NULL,NULL)]
combined.fpkm[,name_or_id:=gene_short_name]
combined.fpkm[is.null(name_or_id),name_or_id:=gene_id]
combined.fpkm[name_or_id=="-" | name_or_id=="",name_or_id:=gene_id]

save(file=output.file,combined.fpkm,star.logs)
