library(data.table)
library(reshape2)
library(ggplot2)
library(grid)

args <- commandArgs(trailingOnly=TRUE)

### ../../data/combined_fpkm
load(args[1])
### ../../data/go_gene_association
load(args[2])
pdf(args[3],width=8,height=8)

hormone.genes <-
    go.gene.association[go_name=="hormone activity" & go_namespace=="molecular_function",gene_name]
data.matrix.dt <- function(x){
    temp <- data.matrix(as.data.frame(x)[,-1])
    rownames(temp) <- as.data.frame(x)[,1]
    temp
}

gene.expression.nozero <- function(genes) {
    temp <-
        data.matrix.dt(dcast(combined.fpkm[human_name %in% genes,],
                             human_name~species,
                             value.var="FPKM",
                             fun.aggregate=function(x){x[1]}))
    temp <-
        temp[apply(temp,1,function(x){sum(x,na.rm=TRUE)})>5,]
    return(temp)
}

pushViewport(viewport(layout=grid.layout(ncol=1,nrow=3)))

### hormone genes
hormone.expression.nozero <-
    gene.expression.nozero(hormone.genes)

hor.exp.long <- melt(log2(hormone.expression.nozero+1)[rev(1:nrow(hormone.expression.nozero)),])
colnames(hor.exp.long)[1:3] <- c("gene","species","fpkm")
print(ggplot(hor.exp.long, aes(y=gene, x=species))
      + geom_tile(aes(fill = fpkm), colour = "white")
      + scale_fill_gradient(low = "white", high = "steelblue")
      + scale_x_discrete("", expand = c(0, 0))
      + scale_y_discrete("", expand = c(0, 0))
      + theme_grey(base_size = 9)
      + theme(legend.position = "none",
              axis.ticks = element_blank(), 
              axis.text.x = element_text(angle = 330, hjust = 0)),
      vp=viewport(layout.pos.col=1,layout.pos.row=1)
      )


### galectin genes
galectin.genes <- unique(grep("^LGALS",combined.fpkm[,human_name],value=TRUE))

galectin.expression.nozero <-
    gene.expression.nozero(galectin.genes)

gal.exp.long <- melt(log2(galectin.expression.nozero+1)[rev(1:nrow(galectin.expression.nozero)),])
colnames(gal.exp.long)[1:3] <- c("gene","species","fpkm")
print(ggplot(gal.exp.long, aes(y=gene, x=species))
      + geom_tile(aes(fill = fpkm), colour = "white")
      + scale_fill_gradient(low = "white", high = "steelblue")
      + scale_x_discrete("", expand = c(0, 0))
      + scale_y_discrete("", expand = c(0, 0))
      + theme_grey(base_size = 9)
      + theme(legend.position = "none",
              axis.ticks = element_blank(), 
              axis.text.x = element_text(angle = 330, hjust = 0)),
      vp=viewport(layout.pos.col=1,layout.pos.row=2)
      )


igfrelated.genes <- unique(grep("^IGF",combined.fpkm[,human_name],value=TRUE))

igfrelated.expression.nozero <-
    gene.expression.nozero(igfrelated.genes)

igf.exp.long <- melt(log2(igfrelated.expression.nozero+1)[rev(1:nrow(igfrelated.expression.nozero)),])
colnames(igf.exp.long)[1:3] <- c("gene","species","fpkm")
print(ggplot(igf.exp.long, aes(y=gene, x=species))
      + geom_tile(aes(fill = fpkm), colour = "white")
      + scale_fill_gradient(low = "white", high = "steelblue")
      + scale_x_discrete("", expand = c(0, 0))
      + scale_y_discrete("", expand = c(0, 0))
      + theme_grey(base_size = 9)
      + theme(legend.position = "none",
              axis.ticks = element_blank(), 
              axis.text.x = element_text(angle = 330, hjust = 0)),
      vp=viewport(layout.pos.col=1,layout.pos.row=3)
      )


popViewport(1)
    
dev.off()
