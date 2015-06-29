args <- c("../../data/combined_fpkm",
          "go_functional_grouping"
          )
args <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(reshape2)
library(ggplot2)
library(grid)


### load("../../data/combined_fpkm")
load(args[1])
### load("go_functional_grouping")
load(args[2])

pdf(file=args[length(args)],width=10,height=20)
data.matrix.dt <- function(x){
    temp <- data.matrix(as.data.frame(x)[,-1])
    rownames(temp) <- as.data.frame(x)[,1]
    temp
}

gene.expression.nozero <- function(genes,expression,min.exp=5,log=FALSE) {
    agg.fun <- function(x){x[1]}
    if (log) {
        agg.fun <- function(x){log2(x[1]+1)}
    }
    temp <-
        data.matrix.dt(dcast(expression[human_name %in% genes,],
                             human_name~species,
                             value.var="FPKM",
                             fun.aggregate=agg.fun))
    return(temp[apply(temp,1,
                      function(x){
                          sum(x,na.rm=TRUE)
                      })> agg.fun(min.exp),])
}



go.terms <- c("^blood vessel development$",
              "^extracellular matrix$",
              "^tube morphogenesis$",
              "apoptosis",
              "lumen$")
go.genes <- list()
for (i  in 1:length(go.terms)) {
    go.term <- go.terms[i]
    go.genes[[go.term]] <-
        go.functional.groups[grepl(go.term,name_1006),hgnc_symbol]
}
# pushViewport(viewport(layout=grid.layout(ncol=1,nrow=length(go.terms),
#                            heights=sapply(go.genes,length))))
for (i  in 1:length(go.terms)) {
    go.term <- go.terms[i]
    go.expression.nozero <-
        gene.expression.nozero(genes=go.genes[[go.term]],
                               expression=combined.fpkm,
                               log=TRUE
                               )
    go.exp.nz.long <-
        melt(go.expression.nozero,
             value.name="fpkm"
             )
    colnames(go.exp.nz.long)[1:2] <- c("gene","species")

    ## reorder by expression level
    go.exp.nz.long$gene <-
        reorder(go.exp.nz.long$gene,
                sapply(go.exp.nz.long$fpkm,function(x){x[is.na(x)] <- 0; x}))

#    pushViewport(viewport(layout.pos.col=1,layout.pos.row=i))
    print(ggplot(go.exp.nz.long,
                 aes(y=gene, x=species))
          + geom_tile(aes(fill = fpkm), colour = "white")
          + scale_fill_gradient()
          + scale_x_discrete("", expand = c(0, 0))
          + scale_y_discrete("", expand = c(0, 0))
          + theme_grey(base_size = 9)
          + theme(legend.position = "none",
                  axis.ticks = element_blank(), 
                  axis.text.x = element_text(angle = 330, hjust = 0))
          + ggtitle(gsub("[\\^\\$]","",go.term))
#          newpage=FALSE
          )
    grid.text(toupper(letters[i]),
              x=unit(0,"npc"),
              y=unit(1,"npc"),
              just="left",
              hjust=1)
#    popViewport()
}
# popViewport()

dev.off()
