
library("data.table")
library("phytools")

### arguments for debugging
args <- c("combined_fpkm","gene_trees","CYP19A1","combined_ancestral_states")

args <- commandArgs(trailingOnly=TRUE)

## combined_fpkm_per_sample
load(args[1])
## gene_trees
load(args[2])

genes <- args[-c(1:2,length(args))]


calculate.ancestral.state <- function(gene,small.distance=0.001) {

    capfirst <- function(s, strict = FALSE) {
        paste(toupper(substring(s, 1, 1)),
        {s <- substring(s, 2); if(strict) tolower(s) else s},
        sep = "")
    }
    ## find the proper tree for this gene
    tree.id <-
        genes.to.tree[combined.fpkm[gene_short_name==gene & species=="homo sapiens",gene_id],
                      unique(tree)]
    if (is.null(tree.id))
        stop(paste0("Unable to find tree for",gene))
    message(paste("gene is",gene,"tree id is",tree.id))
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
    
    ## use the table to label the tree tips
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
    gene.tree.expression[,log_FPKM:=log10(FPKM+1)]
    
    gene.tree.expression <-
        gene.tree.expression[,list(gene_id,gene_short_name,log_FPKM,species)]
    gene.tree.expression[,species:=gsub("^(.)[^ ]* +([^ ][^ ])[^ ]* *$","\\U\\1.\\L\\2",
                                        species,
                                        perl=TRUE
                                        )]
    gene.tree.expression[,duplicate_gene_id:=gene_id]
    
    for (duplicate.gene_id in unique(gene.tree.expression[duplicated(gene_id),gene_id])) {
        num.duplicates <- gene.tree.expression[,sum(gene_id==duplicate.gene_id)]
        duplicate.tree.struct <-
            paste0("(",
                   paste(paste0(rep.int(duplicate.gene_id,num.duplicates),
                                "_",
                                seq.int(1,num.duplicates),
                                ":",small.distance),
                         collapse=","),
                   ")",duplicate.gene_id,":",small.distance,";")
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
            bind.tree(gene.tree.subset,duplicate.tree,where=which(gene.tree.subset$tip==duplicate.gene_id))
    }

    gene.tree.subset <-
        read.tree(text=gsub(":0([\\),])",paste0(":",small.distance,"\\1"),
                            write.tree(gene.tree.subset)))


    gene.tree.expression.vector <-
        gene.tree.expression[,as.numeric(log_FPKM)]
    small.replacement <-
        rexp(n=sum(gene.tree.expression.vector==0),rate=10)
    small.replacement[small.replacement>1e-3] <- 0.1
    gene.tree.expression.vector[gene.tree.expression.vector==0] <-
        small.replacement
    names(gene.tree.expression.vector) <-
        gene.tree.expression[,duplicate_gene_id]
    
    anc.tree <- fastAnc(gene.tree.subset,gene.tree.expression.vector,vars=TRUE)
#    anc.bayes.tree <- anc.Bayes(gene.tree.subset,gene.tree.expression.vector)
    
    return(list(tree=gene.tree.subset,
                expression=gene.tree.expression,
                expression.vector=gene.tree.expression.vector,
                anc.exp=anc.tree$ace,
#                anc.bayes.exp=anc.bayes.tree[nrow(anc.bayes.tree),
#                                             -c(1:2,ncol(anc.bayes.tree))],
                gene=gene)
           )
}

combined.anc.states <-
    lapply(genes,
           calculate.ancestral.state)
names(combined.anc.states) <- genes

save(file=args[length(args)],
     combined.anc.states)
     
