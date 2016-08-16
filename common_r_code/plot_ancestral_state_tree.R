library("data.table")
library("phytools")
library("magrittr")
library("geiger")
##' Plot the Ancestral State Tree
##'
##' See examples in manuscript.Rnw for more details
##' @title plot.ancestral.state.tree
##' @param gene_selector quoted expression to be eval() which selects the genes to include
##' @param genes.to.tree genes.to.tree object
##' @param combined.fpkm combined.fpkm object
##' @param trees trees object 
##' @param small.distance distance to replace 0 with; default is 0.0001
##' @param ... Additional options passed to contMap
##' @return whatever contMap() returns
##' @author Don Armstrong
plot.ancestral.state.tree <- function(gene_selector,genes.to.tree,combined.fpkm,trees,subtree=NULL,small.distance=0.0001) {
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
    if (!is.null(subtree)) {
        gene.tree.mrca <- mrca(gene.tree.subset)
        top.node <- gene.tree.mrca[subtree[1],subtree[2]]
        gene.tree.subset <-
            drop.tip(gene.tree.subset,
                     gene.tree.subset$tip.label[!(gene.tree.subset$tip.label %in%
                                                  tips(gene.tree.subset,
                                                       top.node))])
    }
    gene.tree.expression <-
        gene.tree.expression[gene_id %in% gene.tree.subset$tip.label,]

    gene.tree.expression[,duplicate_gene_id:=gene_id]

    gene.tree.subset <-
        drop.tip(gene.tree.subset,
                 gene.tree.subset$tip.label[!(gene.tree.subset$tip.label %in%
                                              gene.tree.expression[,gene_id])])


    
    for (duplicate.gene_id in unique(gene.tree.expression[duplicated(gene_id),gene_id])) {
        num.duplicates <- gene.tree.expression[,sum(gene_id==duplicate.gene_id)]
        duplicate.tree.struct <-
            paste0("(",
                   paste(paste0(rep.int(duplicate.gene_id,num.duplicates),
                                "_",
                                seq.int(1,num.duplicates),
                                ":",small.distance),
                         collapse=","),
                   "):",small.distance,";")
        ## print(duplicate.tree.struct)
        duplicate.tree <-
            read.tree(text=duplicate.tree.struct)
        ## print(duplicate.tree)
        ## plot(duplicate.tree)
        gene.tree.expression[gene_id==duplicate.gene_id,
                             duplicate_gene_id:=paste0(rep.int(duplicate.gene_id,num.duplicates),
                                                       "_",
                                                       seq.int(1,num.duplicates))]
        gene.tree.subset <-
            bind.tree(gene.tree.subset,
                      duplicate.tree,
                      where=which(gene.tree.subset$tip==duplicate.gene_id))
    }

    ## replace any 0 distances in the tree with small.distance
    gene.tree.subset <-
        read.tree(text=gsub(":0([\\),])",paste0(":",small.distance,"\\1"),
                            write.tree(gene.tree.subset)))

    gene.tree.expression <-
        gene.tree.expression[duplicate_gene_id %in% gene.tree.subset$tip.label,]
    # gene.tree.expression[,log_FPKM:=log10(FPKM+1)]
    
    gene.tree.expression <-
        gene.tree.expression[,list(gene_id,duplicate_gene_id,
                                   gene_short_name,`mean_fpkm`,species)]
    gene.tree.expression[,species:=gsub("^(.)[^ ]* +([^ ][^ ])[^ ]* *$","\\U\\1.\\L\\2",
                                        species,
                                        perl=TRUE
                                        )]
    gene.tree.expression[gene_short_name=="-",
                         gene_short_name:=gsub("ENS[^0-9]+00+0","0..0",
                                               gsub("_\\d+$","",gene_id))]
    gene.tree.expression[,species_gene:=paste0(species," ",gene_short_name)]
    setkey(gene.tree.expression,duplicate_gene_id)
    gene.tree.expression <- gene.tree.expression[gene.tree.subset$tip.label,]
    gene.tree.expression.vector <- gene.tree.expression[,mean_fpkm]
    names(gene.tree.expression.vector) <-
        gene.tree.expression[,species_gene]
    gene.tree.subset$tip.label <- gene.tree.expression[,species_gene]
    theanc <- fastAnc(gene.tree.subset,gene.tree.expression.vector)
    f.tr <- fortify(gene.tree.subset)
    f.tr$fpkm <- NA
    f.tr$fpkm[f.tr$isTip] <-
        as.numeric(gene.tree.expression.vector[f.tr$label[f.tr$isTip]])
    f.tr$fpkm[!f.tr$isTip] <-
        theanc[as.character(f.tr$node[!f.tr$isTip])]
    return(ggtree(f.tr,aes(color=log10(fpkm+1)),size=1)+
           # geom_tiplab(color="black")+
           scale_color_gradient(low="grey90",high="blue",
                                limits=c(0,4),6
                                guide=guide_colorbar(title=expression(log[10](FPKM))))+
           theme(legend.position="bottom"))
}
