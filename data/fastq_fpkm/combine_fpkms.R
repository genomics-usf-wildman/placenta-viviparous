library(data.table)

### this args is for debugging and will be overwritten by the command line arguments
args <- c("protein_id_to_gene_id.txt",
          grep(".fpkm_tracking$",dir(),value=TRUE),
          grep("_trinity_diamond.txt$",dir(),value=TRUE),
          grep("_trinity_align_rsem_isoforms.txt$",dir(),value=TRUE),
          grep("_log_star.txt$",dir(),value=TRUE),
          "combined_fpkms")

args <- commandArgs(trailingOnly=TRUE)

output.file <- args[length(args)]
args <- args[-length(args)]

pb <- txtProgressBar(min=1,max=length(args),style=3)
i <- 0

protein.to.id <- grep("protein_id_to_gene_id.txt",args,value=TRUE)
trinity.rsem <- grep("_trinity_align_rsem_isoforms.txt$",args,value=TRUE)
trinity.diamond.files <- grep("_trinity_diamond.txt$",args,value=TRUE)
gene.files <- grep("genes.fpkm_tracking",args,value=TRUE)
isoform.files <- grep("isoforms.fpkm_tracking",args,value=TRUE)

star.log.files <- grep("_log_star.txt",args,value=TRUE)

## read in the protein to gene id information
protein.to.gene <- fread(protein.to.id)
setkey(protein.to.gene,"protein_id")

## these are the annotations for trinity species
trinity.diamond.species <- NULL
trinity.diamond <- list()
for (file in trinity.diamond.files) {
    species.name <- gsub("_"," ",gsub("_trinity_diamond.txt","",file))
    trinity.diamond.species <-
        c(trinity.diamond.species,
          species.name)
    ## read in appropriate diamond file
    trinity.diamond[[file]] <-
        fread(file)
    setnames(trinity.diamond[[file]],
             c("tracking_id","protein_id","identity.per.len","len",
               "mismatches","gaps","query_begin","query_end",
               "subject_begin","subject_end","evalue","bitscore"))
    ## take the top scoring hit, and use that
    trinity.diamond[[file]] <-
        trinity.diamond[[file]][!duplicated(tracking_id),]
    ## there shouldn't be any bitscores below 40, either.
    trinity.diamond[[file]] <-
        trinity.diamond[[file]][bitscore >= 40,]
    trinity.diamond[[file]][,species:=species.name]
    ## we're going to match on protein_id
    setkey(trinity.diamond[[file]],"protein_id")
    ## map the ensembl ids to gene names
    trinity.diamond[[file]] <-
        protein.to.gene[trinity.diamond[[file]]][,list(tracking_id,gene_id,gene_short_name,
                                                     species)]
    trinity.diamond[[file]][,tracking_id:=paste0(gsub(" ","_",species.name),"_",tracking_id)]
    i <- i+1
    setTxtProgressBar(pb,i)
}

gene.fpkms <- list();
for (file in gene.files) {
    gene.fpkms[[file]] <- fread(file)
    species.name <-
        gsub("_"," ",
             gsub("_genes.fpkm_tracking","",
                  gsub("SRR\\d+_","",file)))
    if (any(species.name %in% trinity.diamond.species)) {
        ## need to merge in the gene names here and then merge down
        ## the calls on the genes to only a single gene
        diamond.file <-
            paste0(gsub(" ","_",species.name),"_trinity_diamond.txt")
        gene.fpkms[[file]][,tracking_id:=paste0(gsub(" +","_",species.name),
                                "_",gsub(":.+","",locus))]
        setkey(gene.fpkms[[file]],"tracking_id")
        setkey(trinity.diamond[[diamond.file]],"tracking_id")
        gene.fpkms[[file]] <-
            trinity.diamond[[diamond.file]][gene.fpkms[[file]]]
        ## sum the FPKM in order to determine the total expression for
        ## this gene
        gene.fpkms[[file]][!is.na(gene_id),
                           FPKM:=sum(FPKM),gene_id]
        gene.fpkms[[file]][!is.na(gene_id),
                           FPKM_conf_lo:=sum(FPKM_conf_lo),gene_id]
        gene.fpkms[[file]][!is.na(gene_id),
                           FPKM_conf_hi:=sum(FPKM_conf_hi),gene_id]
        gene.fpkms[[file]] <-
            gene.fpkms[[file]][is.na(gene_id) | !duplicated(gene_id),]
        gene.fpkms[[file]][,i.gene_id:=NULL]
        gene.fpkms[[file]][,i.gene_short_id:=NULL]
        gene.fpkms[[file]][,i.species:=NULL]
    }
    gene.fpkms[[file]][,species := species.name]
    gene.fpkms[[file]][,file:=file]
    i <- i+1
    setTxtProgressBar(pb,i)
}


gene.fpkms <- rbindlist(gene.fpkms,fill=TRUE)

gene.fpkms[,mean_fpkm := mean(FPKM),tracking_id]
gene.fpkms[,sd_fpkm := sd(FPKM),tracking_id]
gene.fpkms[,n_fpkm := length(FPKM),tracking_id]
# gene.fpkms <- gene.fpkms[!duplicated(tracking_id),]

isoform.fpkms <- list();
for (file in isoform.files) {
    isoform.fpkms[[file]] <- fread(file)
    isoform.fpkms[[file]][,species := gsub("_"," ",
                                       gsub("_isoforms.fpkm_tracking","",
                                            gsub("SRR\\d+_","",file)))]
    isoform.fpkms[[file]][,file:=file]
    i <- i+1
    setTxtProgressBar(pb,i)
}

isoform.fpkms <- rbindlist(c(isoform.fpkms),fill=TRUE)

isoform.fpkms[,mean_fpkm := mean(FPKM,na.rm=TRUE),tracking_id]
isoform.fpkms[,sd_fpkm := sd(FPKM,na.rm=TRUE),tracking_id]
isoform.fpkms[,n_fpkm := length(na.omit(FPKM)),tracking_id]
isoform.fpkms <- isoform.fpkms[!duplicated(tracking_id),]


star.logs <- list()
for (file in star.log.files) {
    star.log <- read.table(file,sep="|",fill=TRUE,stringsAsFactors=FALSE)
    colnames(star.log) <- c("field","value")
    star.log$value <- gsub("\\t","",star.log$value)
    star.log$field <- gsub("(^\\s+|\\s+$)","",star.log$field)
    star.logs[[file]] <- data.table(star.log)[!grepl(":$",field),]
    star.logs[[file]][,species := gsub("_"," ",
                                   gsub("_log_star.txt","",
                                        gsub("SRR\\d+_","",file)))]
    star.logs[[file]][,file:=file]

    i <- i+1
    setTxtProgressBar(pb,i)
}
star.logs <- rbindlist(star.logs)
close(pb)

save(gene.fpkms,
     isoform.fpkms,
     star.logs,
     file=output.file)
