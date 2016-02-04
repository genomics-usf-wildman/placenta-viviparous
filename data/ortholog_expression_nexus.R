library("data.table")
library("reshape2")
library("ape")

args <- c("combined_fpkm","housekeeping_genes_superset",
          "ortholog_expression_matrix.txt")

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
    combined.fpkm.wide[apply(combined.fpkm.wide,1,
                         function(x){sum(as.numeric(x[-1]))>0}),]

combined.fpkm.wide.df <-
    t(data.frame(combined.fpkm.wide[!(human_name %in%
                                      housekeeping.genes.superset[,gene_short_name]),
                                    ][,-1,with=FALSE]))

output.file <- file(args[length(args)],"w")

output <- function(...) {
    cat(file=output.file,
        ...)
}

output("#NEXUS\n")
output("BEGIN DATA;\n")
output(paste0(" dimension ntax=",nrow(combined.fpkm.wide.df),
              " nchar=",ncol(combined.fpkm.wide.df),";\n"))
output("  format datatype=standard symbols=\"0 1\";\n")
output("  matrix\n")

comb.fpkm.joined <- apply(combined.fpkm.wide.df,1,paste,collapse="")

output(paste0(paste("  ",gsub("\\.","_",names(comb.fpkm.joined))," ",comb.fpkm.joined,sep="",collapse="\n"),";\n"))
output("END;\n")

close(output.file)
