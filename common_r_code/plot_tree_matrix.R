plot.tree.matrix <- function(gene.tree,gene.tree.table,expression) {
    gene.tree$tip.label <- gsub("^(ENS.+)\\1$","\\1",gene.tree$tip.label)
    gene.tree.table[,species:=capfirst(gsub("_"," ",species))]
    gene.tree.table[,symbol_or_id:=ifelse(symbol=="NULL",gene_id,paste0(symbol," ",species))]
    gene.tree.expression <-
        combined.fpkm[gene_id %in% gene.tree.table[,gene_id]]
    setkey(gene.tree.table,"protein_id")
    gene.tree$tip.label <-
        gene.tree.table[gene.tree$tip.label,gene_id]
    gene.tree.subset <-
        drop.tip(gene.tree,gene.tree$tip.label[!(gene.tree$tip.label %in%
                                                     gene.tree.expression[,gene_id])])

    gene.tree.expression[,mean_fpkm_log:=log10(mean_fpkm+1)]
    gene.tree.expression.matrix <-
        dcast(gene.tree.expression,
              gene_id~species,value.var="mean_fpkm_log")
    rownames(gene.tree.expression.matrix) <-
        gene.tree.expression.matrix[,1]
    gene.tree.expression.matrix <-
        gene.tree.expression.matrix[,-1]
    colnames(gene.tree.expression.matrix) <-
        capfirst(colnames(gene.tree.expression.matrix))
    colnames(gene.tree.expression.matrix) <-
        gsub("(^)(.)[^ ]*","\\2.",colnames(gene.tree.expression.matrix))
    colnames(gene.tree.expression.matrix) <-
        gsub("( )(..)[^ ]*","\\2.",colnames(gene.tree.expression.matrix))


    setkey(gene.tree.table,"gene_id")
    rownames(gene.tree.expression.matrix) <-
        gene.tree.table[rownames(gene.tree.expression.matrix),symbol_or_id]
        
    
    gene.tree.subset$tip.label <-
        gene.tree.table[gene.tree.subset$tip.label,symbol_or_id]
    p <- ggtree(gene.tree.subset)# %>% add_legend(x=2008, y=5)
    p <- p + geom_tiplab(size=3,align=TRUE)
    
    
    print(gheatmap(p,data=gene.tree.expression.matrix,offset=0.3,colnames=TRUE) +
              theme(legend.title=element_text())
          + scale_fill_gradient(low="white",high="green",na.value="gray50",
                                name=expression(log[10](FPKM)))
          )
}
