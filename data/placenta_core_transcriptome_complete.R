library("data.table")
library("reshape2")
library("org.Hs.eg.db")

args <- c("combined_fpkm","housekeeping_genes_superset","placenta_core_transcriptome_complete")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

combined.fpkm.wide <- 
    data.table(dcast(combined.fpkm,
                     human_name~species,
                     fun.aggregate=sum,
                     value.var="mean_fpkm"))

min.na.zero <- function(x){
    x[is.na(x)] <- 0
    return(min(x))
}

max.na.zero <- function(x){
    return(max(x,na.rm=TRUE))
}

min.but.one <- function(x,n=1){
    x[is.na(x)] <- 0
    while (n > 0) {
        x <- x[-which.min(x)]
        n <- n - 1
    }
    return(min(x))
}

combined.fpkm.wide[,median_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],
                             1,median)]
combined.fpkm.wide[,max_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,max.na.zero)]
combined.fpkm.wide[,all_but_four_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,function(x){min.but.one(x,4)})]
combined.fpkm.wide[,all_but_three_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,function(x){min.but.one(x,3)})]
combined.fpkm.wide[,all_but_two_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,function(x){min.but.one(x,2)})]
combined.fpkm.wide[,all_but_one_expression :=
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,min.but.one)]
combined.fpkm.wide[,all_expression :=
                        apply(combined.fpkm.wide[,-1,with=FALSE],1,min.na.zero)]
combined.fpkm.wide[,all_but_metatherian_expression :=
                        apply(combined.fpkm.wide[,!grepl("^(monodelphis domestica|human_name|.*_expression)$",
                                                         colnames(combined.fpkm.wide)),
                                                 with=FALSE],1,min.na.zero)]
placenta.core.transcriptome.complete <-
    combined.fpkm.wide[,list(human_name,
                             median_expression,
                             all_but_metatherian_expression,
                             all_expression,
                             all_but_one_expression,
                             all_but_two_expression,
                             all_but_three_expression,
                             all_but_four_expression,
                             max_expression)]
setkey(placenta.core.transcriptome.complete,"human_name")
placenta.core.transcriptome.complete[!is.na(human_name),egid:=
  as.vector(unlist(sapply(mget(placenta.core.transcriptome.complete[!is.na(human_name),human_name],
                               org.Hs.egSYMBOL2EG,ifnotfound=NA),
                          function(x){x[1]})))]

save(placenta.core.transcriptome.complete,
     file=args[length(args)])

