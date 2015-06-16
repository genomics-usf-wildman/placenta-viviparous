library(data.table)

args <- commandArgs(trailingOnly=TRUE)

gene.files <- grep("genes.fpkm_tracking",args,value=TRUE)
isoform.files <- grep("isoforms.fpkm_tracking",args,value=TRUE)

gene.fpkms <- list();
for (file in gene.files) {
    gene.fpkms[[file]] <- fread(file)
    gene.fpkms[[file]][,species := gsub("_"," ",
                                    gsub("_genes.fpkm_tracking","",
                                         gsub("SRR\\d+_","",file)))]
    gene.fpkms[[file]][,file:=file]
}

gene.fpkms <- rbindlist(gene.fpkms)

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

isoform.fpkms <- rbindlist(isoform.fpkms)

isoform.fpkms[,mean_fpkm := mean(FPKM),tracking_id]
isoform.fpkms[,sd_fpkm := sd(FPKM),tracking_id]
isoform.fpkms <- isoform.fpkms[!duplicated(tracking_id),]

save(gene.fpkms,isoform.fpkms,
     file=args[length(args)])
