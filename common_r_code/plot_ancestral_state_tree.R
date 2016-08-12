library("data.table")
library("phytools")

plot.ancestral.state.tree <- function(gene_selector,genes.to.tree,combined.fpkm,trees,small.distance=0.0001) {
    capfirst <- function(s, strict = FALSE) {
        paste(toupper(substring(s, 1, 1)),
        {s <- substring(s, 2); if(strict) tolower(s) else s},
        sep = "")
    }
    tree.id <-
        genes.to.tree[combined.fpkm[eval(gene_selector),gene_id],
                      unique(tree)]
    if (length(tree.id)!=1)
        stop(paste0("Number of possible trees is not equal to 1: ",length(tree.id)))
    tree.data <- trees[[tree.id]]$tree.dat
    tree.data <- gsub(":0([\\),])",paste0(":",small.distance,"\\1"),tree.data)
    gene.tree <- read.tree(text=tree.data)
    gene.tree$tip.label <- gsub("^(ENS.+)\\1$","\\1",gene.tree$tip.label)
    
    ## this is the table which corresponds to the gene tree
    gene.tree.table <- trees[[tree.id]]$table
    gene.tree.table[,species:=capfirst(gsub("_"," ",species))]
    gene.tree.table[,short.species:=gsub("^(.)[^ ]* +([^ ]{2})[^ ]*","\\1.\\2.",species)]
    gene.tree.table[,symbol_or_id:=ifelse(symbol=="NULL",gene_id,paste0(symbol," ",short.species))]
    ## use symbol gene_id for duplicates
    gene.tree.table[duplicated(symbol_or_id)|duplicated(symbol_or_id,fromLast=TRUE),
                    symbol_or_id:=paste0(symbol," ",gene_id)]
    setkey(gene.tree.table,"protein_id")
    
    ## convert the protein id into gene_id
    gene.tree$tip.label <-
        gene.tree.table[gene.tree$tip.label,gene_id]

    ## pull expression values for this tree
    gene.tree.expression <-
        combined.fpkm[gene_id %in% gene.tree.table[,gene_id]]
    
    gene.tree.subset <-
        drop.tip(gene.tree,
                 gene.tree$tip.label[!(gene.tree$tip.label %in%
                                       gene.tree.expression[,gene_id])])
    
    gene.tree.expression <-
        gene.tree.expression[gene_id %in% gene.tree.subset$tip.label,]
    # gene.tree.expression[,log_FPKM:=log10(FPKM+1)]
    
    gene.tree.expression <-
        gene.tree.expression[,list(gene_id,gene_short_name,FPKM,species)]
    gene.tree.expression[,species:=gsub("^(.)[^ ]* +([^ ][^ ])[^ ]* *$","\\U\\1.\\L\\2",
                                        species,
                                        perl=TRUE
                                        )]
    gene.tree.expression[gene_short_name=="-",gene_short_name:=gene_id]
    
    contMap(combined.anc.states[[gene]]$tree,
            combined.anc.states[[gene]]$expression.vector)
}
