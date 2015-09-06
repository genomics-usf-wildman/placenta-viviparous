library(data.table)
library(reshape2)
library(ggplot2)
library(grid)
library(scales)

args <- commandArgs(trailingOnly=TRUE)

### ../../data/combined_fpkm
load(args[1])
### ../../data/go_gene_association
load(args[2])
### ../species_ordering.R
source(args[3])


plot_type <- "onepage"
if (any(grepl("multipage",args[length(args)]))) {
    plot_type <- "multipage"
}
pdf(args[length(args)],width=8,height=8)

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

### hormone genes
hormone.genes <- c("CGA","CGB","CGB1","CGB2","CGB5","CGB8","CSH1","CSH2","CSHL1","CYP19A1",
                   "GH2","GHRH","HSD3B1","NPPB","NPPC","PRL","RLN1")
hormone.expression.nozero <-
    gene.expression.nozero(hormone.genes)



hor.exp.long <- melt(log2(hormone.expression.nozero+1)[rev(1:nrow(hormone.expression.nozero)),])
colnames(hor.exp.long)[1:3] <- c("gene","species","fpkm")
hor.exp.long <- hor.exp.long[hor.exp.long$gene!="IGF2",]

### galectin genes
galectin.genes <- unique(grep("^LGALS",combined.fpkm[,human_name],value=TRUE))

galectin.expression.nozero <-
    gene.expression.nozero(galectin.genes)

gal.exp.long <- melt(log2(galectin.expression.nozero+1)[rev(1:nrow(galectin.expression.nozero)),])
colnames(gal.exp.long)[1:3] <- c("gene","species","fpkm")

### IGF related genes

igfrelated.genes <- unique(grep("^IGF",combined.fpkm[,human_name],value=TRUE))

igfrelated.expression.nozero <-
    gene.expression.nozero(igfrelated.genes)

igf.exp.long <- melt(log2(igfrelated.expression.nozero+1)[rev(1:nrow(igfrelated.expression.nozero)),])
colnames(igf.exp.long)[1:3] <- c("gene","species","fpkm")

scalerange <- range(c(igf.exp.long$fpkm,gal.exp.long$fpkm,hor.exp.long$fpkm),na.rm=TRUE)
gradientends <- scalerange + rep(c(0,100,200), each=2)
colorends <- c("white", "red", "white", "green", "white", "blue")

gal.exp.long$fpkm <- gal.exp.long$fpkm+100
hor.exp.long$fpkm <- hor.exp.long$fpkm+200

combined.long <-
    rbind(igf.exp.long,
          gal.exp.long,
          hor.exp.long)
combined.long$species <-
    factor(combined.long$species,
           levels=species.ordering)

common.plot.options <- 
    list(geom_tile(aes(fill = fpkm), colour = "white"),
         scale_fill_gradientn(colours = colorends, values = rescale(gradientends)),
         ##   + scale_fill_gradient(low = "white", high = "steelblue")
         scale_x_discrete("", expand = c(0, 0)),
         scale_y_discrete("", expand = c(0, 0)),
         theme_grey(base_size = 9),
         theme(legend.position = "none",
               axis.ticks = element_blank(), 
               axis.text.x = element_text(angle = 300, hjust = 0, vjust=1))
         )
if (plot_type=="multipage") {
    print(ggplot(combined.long[1:nrow(igf.exp.long),],aes(y=gene,x=species))
          + common.plot.options
          + ggtitle("IGF Related Genes")
          )
    print(ggplot(combined.long[nrow(igf.exp.long)+1:nrow(gal.exp.long),],aes(y=gene,x=species))
          + common.plot.options
          + ggtitle("Galectins")
          )
    print(ggplot(combined.long[nrow(igf.exp.long)+nrow(gal.exp.long)+1:nrow(hor.exp.log),],
                 aes(y=gene,x=species))
          + common.plot.options
          + ggtitle("Hormones")
          )
} else {
    print(ggplot(combined.long, aes(y=gene, x=species))
          + common.plot.options
          )
}

dev.off()
