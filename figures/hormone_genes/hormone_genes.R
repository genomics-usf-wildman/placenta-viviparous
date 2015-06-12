library(data.table)
library(reshape2)
library(ggplot2)
library(heatmap3)

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
pdf(args[3],width=8,height=8)

hormone.genes <-
    go.gene.association[go_name=="hormone activity" & go_namespace=="molecular_function",gene_name]
data.matrix.dt <- function(x){
    temp <- data.matrix(as.data.frame(x)[,-1])
    rownames(temp) <- as.data.frame(x)[,1]
    temp
}
hormone.expression <-
    data.matrix.dt(dcast(combined.fpkm[human_name %in% hormone.genes,],
                         human_name~species,
                         value.var="FPKM",
                         fun.aggregate=function(x){x[1]}))

hormone.expression.nozero <-
    hormone.expression[apply(hormone.expression,1,function(x){sum(x,na.rm=TRUE)})>5,]
hormone.expression.nozero[is.na(hormone.expression.nozero)] <- 0

heatmap3(data.matrix(hormone.expression.nozero),
         scale=NULL)
dev.off()
