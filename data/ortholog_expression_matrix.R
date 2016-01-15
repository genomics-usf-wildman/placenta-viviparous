library("data.table")
library("reshape2")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

combined.fpkm.wide <- 
    data.table(dcast(combined.fpkm,
                     human_name~species,
                     fun.aggregate=function(x){as.numeric(sum(x)>1)},
                     value.var="mean_fpkm"))

combined.fpkm.wide <-
    combined.fpkm.wide[!is.na(human_name),]

combined.fpkm.wide <-
    combined.fpkm.wide[!(human_name %in% housekeeping.genes.superset[,gene_short_name]),]


write.table(file=args[length(args)],
            sep="\t",
            combined.fpkm.wide)
