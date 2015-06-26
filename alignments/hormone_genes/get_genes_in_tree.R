library(data.table)
args <- commandArgs(trailingOnly=TRUE)
load(args[3])
load(args[4])

gene.name <- args[1]
setkey(genes.to.tree,"gene_id")

get.genes.in.tree <- function(gene.name) {
    gene.id <- combined.fpkm[gene_short_name==gene.name &
                                 species=="homo sapiens",gene_id]
    tree.id <- genes.to.tree[gene.id,tree]
    return(trees[[tree.id]]$table[,gene_id])
}

genes <-
    intersect(get.genes.in.tree(gene.name),
              combined.fpkm[,gene_id])

cat(genes,sep="\n",file=args[2])

