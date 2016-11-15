
library(reshape2)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

### combined_fpkm
load(args[1])

homology.table <-
    dcast(combined.fpkm,
          human_name~species,
          value.var="mean_fpkm",
          fun.aggregate=function(x){x[1]})

### ditch rows which are all zero

gene.sum <-
    apply(homology.table,1,function(x){sum(as.numeric(x[-1]),na.rm=TRUE)})
homology.table <- homology.table[gene.sum>0 & !is.na(homology.table$human_name),]
rownames.homology.table <-
    homology.table[,1]
homology.table <- data.matrix(homology.table[,-1])
rownames(homology.table) <-
    rownames.homology.table

save(homology.table,
     file=args[length(args)])
