library(data.table)
library(reshape2)
library(ggplot2)
library(heatmap3)

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
pdf(args[3],width=8,height=8)

combined.fpkm[,top.gene:=FALSE]

### lets ignore mt genes
combined.fpkm <-
    combined.fpkm[is.na(gene_short_name) |
                      !(grepl("^MT-",gene_short_name)|grepl("^MT-",human_name)),]
combined.fpkm <-
    combined.fpkm[is.na(gene_short_name) |
                      !(grepl("^mt-",gene_short_name)|grepl("^mt-",human_name)),]
combined.fpkm <-
    combined.fpkm[is.na(gene_short_name) |
                      !(grepl("rRNA",gene_short_name)|grepl("rRNA",human_name)),]

for (specie in levels(factor(combined.fpkm[,species]))) {
    top.genes <- combined.fpkm[species==specie,top.gene]
    top.genes[order(combined.fpkm[species==specie,mean_fpkm],decreasing=TRUE)[1:20]] <- TRUE
    combined.fpkm[species==specie,top.gene:=top.genes]
}
combined.fpkm[!is.na(human_name) &
                  human_name %in% combined.fpkm[top.gene==TRUE,human_name],
              top.gene:=TRUE]

top.genes <-
    combined.fpkm[combined.fpkm[,top.gene],]

top.genes[,human_name_or_species:=human_name]

top.genes[is.na(human_name_or_species) &
              !is.na(gene_short_name) &
                  !gene_short_name == "-",
          human_name_or_species:=paste(sep="-",species,gene_short_name)]

top.genes[is.na(human_name_or_species),
          human_name_or_species:=tracking_id]

data.matrix.dt <- function(x){
    temp <- data.matrix(as.data.frame(x)[,-1])
    rownames(temp) <- as.data.frame(x)[,1]
    temp
}
top.genes.table <-
    data.matrix.dt(dcast(top.genes,
                         human_name_or_species~species,
                         value.var="mean_fpkm",
                         fun.aggregate=function(x){x[1]}))

### order the columns by the dendrogram
top.genes.table <-
    top.genes.table[,c(1,7,3,8,4,6,5,2)]
### reorder top.genes.table by species in which it is top
top.genes.table <-
    top.genes.table[order(unname(apply(top.genes.table,1,which.max)),
                          unname(apply(top.genes.table,1,function(x){max(x,na.rm=TRUE)})),
                          decreasing=TRUE),]


heatmap3(log2(t(top.genes.table)+1),Colv=NA,
         distfun=dist,
         scale="none",
         margins=c(5,7),
         ylab="Species",
         xlab="Genes",
         main="Top 20 Genes in Each Species",
         cexCol=0.35
         )

dev.off()
