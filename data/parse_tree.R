library(ape)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)


parse.tree <- function(texttree) {
    ## remove extended newick since ape's read tree doesn't like it
    the.tree <- read.tree(text=gsub("\\[[^\\]]+\\]","",texttree))
    the.tree$tip.label <-
        gsub("^(.{5,})\\1$","\\1",the.tree$tip.label)
    return(the.tree)
}


file.con <- file(args[1])
open(file.con)

trees <- list()
## because the tree ID isn't present in the file (you have to search
## the root tree db which is kind of lame), we'll just use numbers for
## now.
current.tree <- 1
next.data <- FALSE
read.size <- 1000000
while (1) {
    results <- readLines(file.con,n=read.size,ok=TRUE,warn=FALSE);
    for (i in 1:length(results)) {
        if(length(trees) < current.tree) {
            trees[[current.tree]] <-
                list(seqs=vector(mode="character"),
                     tree.dat=NULL,
                     tree=NULL)
        }
        if (next.data) {
            trees[[current.tree]]$tree.dat <-
                results[i]
                ## parse.tree(results[i])
            next.data <- FALSE
            next;
        }
        if (any(grepl("^SEQ",results[i]))) {
            trees[[current.tree]]$seqs[length(trees[[current.tree]]$seqs)+1] <-
                gsub("^SEQ\ +","",results[i])
            next;
        }
        if (any(grepl("^DATA",results[i]))) {
            next.data <- TRUE
            next;
        }
        if (any(grepl("^//",results[i]))) {
            current.tree <- current.tree + 1
            next;
        }
    }
    print(current.tree)
    print(length(results))
    print(seek(file.con))
    if (length(results) < read.size) {
        break;
    }
}
close(file.con)

genes.to.tree <- list()
prot.to.tree <- list()
for (current.tree in 1:length(trees)) {
    if (length(trees[[current.tree]]$seq) < 1)
        next;
    trees[[current.tree]]$table <-
        fread(paste(sep="",collapse="\n",
                    gsub(" \\(\\d+ of \\d+\\)$","",
                         trees[[current.tree]]$seq),""))
    setnames(trees[[current.tree]]$table,
             c("species","protein_id","chr","start","stop","strand","gene_id","symbol"))
    genes.to.tree[[current.tree]] <-
        data.table(gene_id=trees[[current.tree]]$table[,gene_id],
                   tree=current.tree)
    prot.to.tree[[current.tree]] <-
        data.table(protein_id=trees[[current.tree]]$table[,protein_id],
                   tree=current.tree)
}

genes.to.tree <-
    rbindlist(genes.to.tree)
setkey(genes.to.tree,"gene_id")
prot.to.tree <-
    rbindlist(prot.to.tree)
setkey(prot.to.tree,"protein_id")

save(file=args[length(args)],
     trees,
     genes.to.tree,
     prot.to.tree)
