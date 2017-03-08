library("data.table")
library("reshape2")
library("parallel")

args <- commandArgs(trailingOnly=TRUE)

if (is.null(getOption("mc.cores"))) {
    options(mc.cores=12)
}

output.file <- args[1]

## fastq_fpkm/combined_fpkms
load(args[3])

## oma-groups_long.txt
oma.groups <- fread(args[2])
oma.groups[,species:=gsub("\\d+","",oma_entry)]
oma.groups <- oma.groups[grepl("G0",ensembl),]
setkey(oma.groups,"ensembl")
oma.groups.human <- oma.groups[species=="HUMAN"]
setkey(oma.groups.human,"fingerprint")


### ensembl_all_orthologs
load(args[4])
orthologs.subset <-
    orthologs[secondary=="hsapiens",]
setkey(orthologs.subset,"ensembl_gene_id")

human_gene_names <-
    gene.fpkms[species=="homo sapiens",list(gene_short_name=gene_short_name[1]),by=gene_id]
setkey(human_gene_names,"gene_id")

identify.homolog <- function(ensembl.id) {
    ## identify the group number from OMA
    oma.fingerprint <- oma.groups[ensembl.id,fingerprint]
    ## if there are any human genes with this fingerprint, great, use
    ## it
    if (!is.null(oma.fingerprint) &
        !is.na(oma.fingerprint)) {
        human.id <- oma.groups.human[oma.fingerprint,
                                     ensembl]
        if (!is.null(human.id) && !is.na(human.id)) {
            homolog <-
                human_gene_names[human.id,gene_short_name][1]
            if(!is.null(homolog)) {
                return(homolog)
            }
        }
    }
    ## if not, if there is no ensembl ortholog, return NA
    hum.id <- orthologs.subset[ensembl.id,]
    if (is.null(hum.id) || is.na(hum.id[,homolog_ensembl_gene_id])) {
        return(NA)
    }
    ## if that ortholog is a 1:1 ortholog, and there isn't another OMA
    ## group which is that ortholog, use it, otherwise it's the wrong
    ## ortholog
    if (hum.id[,all(homolog_orthology_type=="ortholog_one2one")]) {
        oma.species <- oma.groups[ensembl.id,species]
        hum.gene.id <- hum.id[,homolog_ensembl_gene_id]
        human.oma.fingerprint <- oma.groups[hum.gene.id,fingerprint]
        if (is.na(human.oma.fingerprint) ||
            length(oma.groups[species==oma.species & fingerprint==human.oma.fingerprint,ensembl])==0 ||
            is.na(oma.groups[species==oma.species & fingerprint==human.oma.fingerprint,ensembl])) {
            return(human_gene_names[hum.gene.id,
                                    gene_short_name[1]])
        } else {
            return(NA)
        }
    }
    ## if there is a 1:many ortholog, things are more complicated. If
    ## there is an existing OMA group with a human member, then return NA.
    hum.many <-
        orthologs.subset[homolog_ensembl_gene_id %in% hum.id[,homolog_ensembl_gene_id] &
                         primary %in% hum.id[,primary],]
    this.oma.group <-
        oma.groups[hum.many[,ensembl_gene_id],]
    if (!all(is.null(oma.groups.human[this.oma.group[!is.na(fingerprint),fingerprint],ensembl])) &&
        !all(is.na(oma.groups.human[this.oma.group[!is.na(fingerprint),fingerprint],ensembl]))) {
        return(NA)
    }
    ## if no existing OMA group has a human member, but there is only
    ## one OMA group, that 1:many ortholog is the right one
    if (this.oma.group[,sum(!is.na(fingerprint))==1] &&
        this.oma.group[,any(ensembl==ensembl.id & !is.na(fingerprint))]) {
        hum.gene.id <- hum.id[,homolog_ensembl_gene_id]
        return(human_gene_names[hum.gene.id,
                                gene_short_name[1]])
    }
    ## otherwise, this is a 1:many ortholog for which we don't know enough about; return NA
    return(NA)
}

### double check that the homolog function is working properly
if (!(identify.homolog("ENSDNOG00000004237")=="ANXA2" &&
      is.na(identify.homolog("ENSDNOG00000048444")) &&
      is.na(identify.homolog("ENSBTAG00000000025")) &&
      identify.homolog("ENSSSCG00000022230")=="CD9" &&
      is.na(identify.homolog("ENSSSCG00000000710")))
    ) {
    print(identify.homolog("ENSDNOG00000004237"))
    stop("Something is wrong with identify.homolog")
}

gene.fpkms[species=="homo sapiens",human_name:=gene_short_name]
gene.fpkms[human_name=="",human_name:=NA]
gene.fpkms[human_name=="-",human_name:=NA]
gene.fpkms[grepl(" ",human_name),
              human_name:=NA]

setkey(gene.fpkms,"tracking_id")
gene.fpkms[is.na(human_name) & !is.na(gene_id),
           human_name := mcmapply(identify.homolog,gene_id)]

## remove all non-1:1 orthologs
non.1.1.orthologs <-
     gene.fpkms[,list(gene_id,n.orthologs=length(gene_id)),
                by=list(human_name,file)][n.orthologs>1,unique(gene_id)]
if (!is.null(non.1.1.orthologs)) {
    gene.fpkms[non.1.1.orthologs,
               human_name:=NA]
}

### this shouldn't be required, but I've triggered some bug in the
### data.table code
combined.fpkm <- data.table(data.frame(gene.fpkms))

setkey(oma.groups,"ensembl")
setkey(combined.fpkm,"tracking_id")
combined.fpkm <- oma.groups[,list(group_num,fingerprint,ensembl)][combined.fpkm]
setnames(combined.fpkm,"ensembl","tracking_id")
setnames(combined.fpkm,"group_num","oma_group_num")

combined.fpkm[,oma_group_name:=as.character(NA)]
combined.fpkm[!is.na(oma_group_num),
              oma_group_name:=paste(collapse=",",
                  grep("^-$",unique(gene_short_name),value=TRUE,invert=TRUE)),
              by=as.numeric(oma_group_num)]

combined.fpkm[,c("class_code","nearest_ref_id",
                 "i.gene_short_name",
                 "tss_id","locus","length","coverage"):=
                     list(NULL,NULL,NULL,NULL,NULL,NULL,NULL)]
combined.fpkm[,name_or_id:=gene_short_name]
combined.fpkm[is.null(name_or_id),name_or_id:=gene_id]
combined.fpkm[name_or_id=="-" | name_or_id=="",name_or_id:=gene_id]

### this is wrong for pan paniscus, so fixing it here
combined.fpkm[human_name=="GH1" & oma_group_name=="GH2",human_name:=NA]

### diamond gene_to_human_gene.txt
diamond.align <- fread(args[5])
setkey(diamond.align,"gene_id")
setnames(diamond.align,"human_gene_symbol","human_alignment_symbol")
setkey(combined.fpkm,"gene_id")
combined.fpkm <-
    diamond.align[,list(gene_id,human_alignment_symbol)][combined.fpkm]

### gtf types
gtf.types <- fread(args[6])
setkey(gtf.types,"gene_id")
setkey(combined.fpkm,"gene_id")
combined.fpkm <-
    gtf.types[combined.fpkm]
setkey(combined.fpkm,"tracking_id")

combined.fpkm[species=="homo sapiens", human_alignment_symbol:=human_name]
combined.fpkm[species=="ateles fusciceps", human_alignment_symbol:=human_name]

if (!all(grepl("_per_sample",output.file))) {
    combined.fpkm <-
        combined.fpkm[!duplicated(tracking_id)]
}

save(file=output.file,combined.fpkm,star.logs)
