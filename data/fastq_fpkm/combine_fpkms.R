library(data.table)

### this args is for debugging and will be overwritten by the command line arguments
args <- c("protein_id_to_gene_id.txt",
          grep(".fpkm_tracking$",dir(),value=TRUE),
          grep("_trinity_diamond.txt$",dir(),value=TRUE),
          grep("_trinity_align_rsem_isoforms.txt$",dir(),value=TRUE),
          "combined_fpkms")

args <- commandArgs(trailingOnly=TRUE)

output.file <- args[length(args)]
args <- args[-length(args)]

protein.to.id <- grep("protein_id_to_gene_id.txt",args,value=TRUE)
trinity.rsem <- grep("_trinity_align_rsem_isoforms.txt$",args,value=TRUE)
trinity.diamond.files <- grep("_trinity_diamond.txt$",args,value=TRUE)
gene.files <- grep("genes.fpkm_tracking",args,value=TRUE)
isoform.files <- grep("isoforms.fpkm_tracking",args,value=TRUE)

## read in the protein to gene id information
protein.to.gene <- fread(protein.to.id)
setkey(protein.to.gene,"protein_id")

trinity.fpkms <- list()
trinity.gene.fpkms <- list()
## these are the annotations
trinity.diamond <- list()
for (file in trinity.rsem) {
    species.name <- gsub("_trinity_align_rsem_isoforms.txt","",file)
    diamond.file <- paste0(species.name,"_trinity_diamond.txt")
    
    if (!any(diamond.file %in% trinity.diamond.files)) {
        stop(paste0("There is no diamond file corresponding to species ",species))
    }
    trinity.fpkms[[file]] <- fread(file)
    trinity.fpkms[[file]][,species:=gsub("_"," ",species.name)]
    setnames(trinity.fpkms[[file]],"transcript_id","tracking_id")
    ## read in appropriate diamond file
    trinity.diamond[[file]] <-
        fread(diamond.file)
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
    setkey(trinity.diamond[[file]],"tracking_id")
    setkey(trinity.fpkms[[file]],"tracking_id")
    trinity.fpkms[[file]] <-
        trinity.diamond[[file]][trinity.fpkms[[file]]]
    setkey(trinity.fpkms[[file]],"protein_id")
    protein.to.gene[trinity.fpkms[[file]]]
    trinity.fpkms[[file]] <-
        protein.to.gene[trinity.fpkms[[file]]][,list(tracking_id,gene_id,gene_short_name,FPKM)]
    trinity.fpkms[[file]][,tracking_id:=paste0(species.name,"_",tracking_id)]
    trinity.fpkms[[file]][,file:=file]
    trinity.gene.fpkms[[file]] <-
        copy(trinity.fpkms[[file]])
    ## I think summing the FPKM per gene is the best approach here
    trinity.gene.fpkms[[file]][,FPKM.sum:=sum(FPKM),gene_id]
    trinity.gene.fpkms[[file]] <-
        trinity.gene.fpkms[[file]][!is.na(gene_id) & !duplicated(gene_id),]
}

gene.fpkms <- list();
for (file in gene.files) {
    gene.fpkms[[file]] <- fread(file)
    gene.fpkms[[file]][,species := gsub("_"," ",
                                    gsub("_genes.fpkm_tracking","",
                                         gsub("SRR\\d+_","",file)))]
    gene.fpkms[[file]][,file:=file]
}

gene.fpkms <- rbindlist(c(gene.fpkms,trinity.gene.fpkms),fill=TRUE)

gene.fpkms[,mean_fpkm := mean(FPKM),tracking_id]
gene.fpkms[,sd_fpkm := sd(FPKM),tracking_id]
gene.fpkms <- gene.fpkms[!duplicated(tracking_id),]

isoform.fpkms <- list();
for (file in isoform.files) {
    isoform.fpkms[[file]] <- fread(file)
    isoform.fpkms[[file]][,species := gsub("_"," ",
                                       gsub("_isoforms.fpkm_tracking","",
                                            gsub("SRR\\d+_","",file)))]
    isoform.fpkms[[file]][,file:=file]
}

isoform.fpkms <- rbindlist(c(isoform.fpkms,trinity.fpkms),fill=TRUE)

isoform.fpkms[,mean_fpkm := mean(FPKM),tracking_id]
isoform.fpkms[,sd_fpkm := sd(FPKM),tracking_id]
isoform.fpkms <- isoform.fpkms[!duplicated(tracking_id),]

save(gene.fpkms,
     isoform.fpkms,
     file=output.file)
