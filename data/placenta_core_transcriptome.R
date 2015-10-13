library("data.table")
library("reshape2")
library("org.Hs.eg.db")

args <- c("combined_fpkm","housekeeping_genes_superset","placenta_core_transcriptome")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

combined.fpkm.wide <- 
    data.table(dcast(combined.fpkm,
                     human_name~species,
                     fun.aggregate=sum,
                     value.var="mean_fpkm"))

combined.fpkm.wide[,median_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],
                             1,median)]
combined.fpkm.wide[,min_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,min)]
setkey(combined.fpkm.wide,"human_name")

## all genes which are not housekeeping genes which have minimum 
core.placenta.transcriptome.long <- 
    data.table(melt(combined.fpkm.wide[min_expression >= 10,
                                       ][order(-median_expression)
                                         ][!(human_name %in% 
                                                 housekeeping.genes.superset[,gene_short_name])],
                    id.vars="human_name",
                    variable.name="Species",
                    value.name="FPKM"
                    ))
core.placenta.transcriptome.long[,human_name:=
                                     reorder(human_name,
                                             combined.fpkm.wide[core.placenta.transcriptome.long[,human_name],median_expression])]


all.genes <- combined.fpkm.wide[!is.na(human_name),unique(human_name)]
placenta.transcriptome.genes <- 
    core.placenta.transcriptome.long[!is.na(human_name),unique(human_name)]


pts <- data.table(all.genes=all.genes)
pts <- pts[,placenta.ts:=all.genes %in% placenta.transcriptome.genes]
pts <- pts[,housekeeping:=all.genes %in% housekeeping.genes.superset[,gene_short_name]]
pts <- pts[,egid:=as.vector(unlist(mget(all.genes,org.Hs.egSYMBOL2EG,ifnotfound=NA)))]
## this list contains all possible genes with a 0 if it is not in the
## core placenta transcriptome, or a 1 if it is in the core placenta
## transcriptome
pts.list <- factor(as.integer(pts[!is.na(egid) & ! housekeeping,placenta.ts]))
names(pts.list) <- pts[!is.na(egid) & ! housekeeping,egid]

placenta.transcriptome.list <-
    pts.list

save(core.placenta.transcriptome.long,
     placenta.transcriptome.list,
     pts,
     file=args[length(args)])


