library("magrittr")
library("tidyr")
library("geiger")


don.gheatmap <- function(p, data, offset=0, width=1, low="green", high="red",
                     color="white", colnames=TRUE, colnames_position="bottom", font.size=4,font.angle=90) {

    colnames_position %<>% match.arg(c("bottom", "top"))
    variable <- value <- lab <- y <- NULL
    
    ## if (is.null(width)) {
    ##     width <- (p$data$x %>% range %>% diff)/30
    ## }

    ## convert width to width of each cell
    width <- width * (p$data$x %>% range %>% diff) / ncol(data)
    
    isTip <- x <- y <- variable <- value <- from <- to <- NULL
 
    df=p$data
    df=df[df$isTip,]
    start <- max(df$x) + offset

    dd <- data[df$label[order(df$y)],]
    dd$y <- sort(df$y)

    dd$lab <- rownames(dd)
    ## dd <- melt(dd, id=c("lab", "y"))
    dd <- gather(dd, variable, value, -c(lab, y))
    
    dd$value[is.na(dd$value) | dd$value == ""] <- NA

    V2 <- start + as.numeric(dd$variable) * width
    mapping <- data.frame(from=dd$variable, to=V2)
    mapping <- unique(mapping)

    dd$x <- V2

    #if (is.null(color)) {
        p2 <- p + geom_raster(data=dd, aes(x, y, fill=value), inherit.aes=FALSE)
##    } else {
##        p2 <- p + geom_raster(data=dd, aes(x, y, fill=value), color=color, inherit.aes=FALSE)
##    }
    if (is(dd$value,"numeric")) {
        p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value="white")
    } else {
        p2 <- p2 + scale_fill_discrete(na.value="white")
    }
    
    if (colnames) {
        if (colnames_position == "bottom") {
            y <- min(p$data$y)-0.6
        } else {
            y <- max(p$data$y) + 1
        }
        p2 <- p2 + geom_text(data=mapping,
                             aes(x=to, label=from),
                             angle=font.angle,
                             hjust=ifelse(font.angle==0,0.5,1),
                             vjust=0.5,
                             y=y, size=font.size, fontface="italic",inherit.aes = FALSE)
    }

    p2 <- p2 + theme(legend.position="top",
                     legend.title=element_blank(),
                     legend.margin=unit(0,"mm"),
                     plot.margin=unit(c(0,0,0,0),"mm")
                     )
    # p2 <- p2 + guides(fill = guide_legend(override.aes = list(color = NULL)))
    
    attr(p2, "mapping") <- mapping
    return(p2)
}

plot.tree.matrix <- function(gene.tree,gene.tree.table,gene.tree.expression,offset.ratio=0.5,min.fpkm=NA,subtree=NULL,vp=NULL,fontsize=3,width=1.5,axis.text.size=5,axis.text.angle=90) {
    gene.tree$tip.label <- gsub("^(ENS.+)\\1$","\\1",gene.tree$tip.label)
    gene.tree.table[,species:=capfirst(gsub("_"," ",species))]
    gene.tree.table[,short.species:=gsub("^(.)[^ ]* +([^ ]{2})[^ ]*","\\1.\\2.",species)]
    gene.tree.table[,symbol_or_id:=ifelse(symbol=="NULL",gene_id,paste0(symbol," ",short.species))]
    ## use symbol gene_id for duplicates
    gene.tree.table[duplicated(symbol_or_id)|duplicated(symbol_or_id,fromLast=TRUE),
                    symbol_or_id:=paste0(symbol," ",gene_id)]
    gene.tree.expression.subset <- gene.tree.expression
    if (!is.na(min.fpkm)) {
        gene.tree.expression.subset <-
            gene.tree.expression[mean_fpkm >= min.fpkm ,]
    }
    setkey(gene.tree.table,"protein_id")
    gene.tree$tip.label <-
        gene.tree.table[gene.tree$tip.label,gene_id]
    gene.tree.subset <-
        drop.tip(gene.tree,
                 gene.tree$tip.label[!(gene.tree$tip.label %in%
                                       gene.tree.expression.subset[,gene_id])])
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
        gsub("( +)([^ ][^ ])[^ ]* *$","\\2",colnames(gene.tree.expression.matrix))


    setkey(gene.tree.table,"gene_id")
    rownames(gene.tree.expression.matrix) <-
        gene.tree.table[rownames(gene.tree.expression.matrix),symbol_or_id]
        
    
    gene.tree.subset$tip.label <-
        gene.tree.table[gene.tree.subset$tip.label,symbol_or_id]
    p <- ggtree(gene.tree.subset)# %>% add_legend(x=2008, y=5)
#    p <- p; + geom_tiplab(aes(label=""),align=TRUE,size=10)
    

    text.labels <- p$data[!is.na(p$data$label),]
    text.labels$xpos <- max(p$data$x)+max(p$data$x)*offset.ratio
    text.labels$face <- "italic"
    text.labels$face[grepl("^ENS",text.labels$label)] <- "plain"
    g1 <-
        (don.gheatmap(p,data=gene.tree.expression.matrix,
                      width=width,
                      offset=max(p$data$x)*offset.ratio,
                      font.size=axis.text.size,
                      font.angle=axis.text.angle,
                      colnames=TRUE)
            + geom_text(aes(x=xpos,
                            y=y,
                            label=label,
                            hjust=1,
                            fontface=face
                            ),
                        size=fontsize,
                        data=text.labels
                        )
            + theme(legend.title=element_text(),legend.position="top",
                    legend.margin=unit(0,"mm"),
                    plot.margin=unit(c(2,0,0,0),"mm"),
                    panel.grid=element_line(size=0.5,color="gray50"),
                    axis.text=element_text(size=axis.text.size,face="italic")
                    )
            + scale_fill_gradient(low = "#132B43", high = "#56B1F7",
                                  space = "Lab", na.value = NA,
                                  name=expression(log[10](FPKM)))
        )
    g1 <- ggplot_gtable(ggplot_build(g1))
    g1$layout$clip[g1$layout$name=="panel"] <- "off"
    g1$vp=vp
    grid.draw(g1)
}

