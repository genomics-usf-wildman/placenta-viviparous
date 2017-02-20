library("data.table")
library("reshape2")
library("org.Hs.eg.db")
library("Biostrings")

args <- c("combined_fpkm","gene_trees","housekeeping_genes_superset","placenta_core_transcriptome")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
load(args[3])

combined.fpkm <- genes.to.tree[combined.fpkm]

collapse.gene.names <- function(x,min.collapse=2) {
    ## longest common substring
    if (is.null(x) || length(x)==0) {
        return(as.character(NA))
    }
    x <- sort(unique(x))
    str_collapse <- function(y,len) {
        if (len == 1 || length(y) < 2) {
            return(y)
        }
        y.tree <-
            gsub(paste0("^(.{",len,"}).*$"),"\\1",y[1])
        y.rem <-
            gsub(paste0("^.{",len,"}"),"",y)
        y.rem.prefix <-
            sum(combn(y.rem,2,function(x){lcprefix(x[1],x[2])}) >= 2)
        if (length(y.rem) > 3 &&
            y.rem.prefix >= 2
            ) {
            y.rem <- 
                collapse.gene.names(y.rem,min.collapse=1)
        }
        paste0(y.tree,
               "{",paste(collapse=",",
                         y.rem),"}")
    }
    i <- 1
    ret <- NULL
    while (i <= length(x)) {
        col.pmin <-
            pmin(sapply(x,lcprefix,x[i]))
        collapseable <-
            which(col.pmin > min.collapse)
        if (length(collapseable) == 0) {
            ret <- c(ret,x[i])
            i <- i+1
        } else {
            ret <- c(ret,
                     str_collapse(x[collapseable],
                                  min(col.pmin[collapseable]))
                     )
            i <- max(collapseable)+1
        }
    }
    return(paste0(collapse=",",ret))
}

combined.fpkm[!is.na(tree),
              tree_genes:=collapse.gene.names(na.omit(human_name)),
              by=tree]

combined.fpkm.wide <-
    data.table(dcast(combined.fpkm[!is.na(tree),],
                     tree+tree_genes~species,
                     fun.aggregate=mean,
                     value.var="mean_fpkm"))

min.na.zero <- function(x){
    x[is.na(x)] <- 0
    return(min(x))
}

max.na.zero <- function(x){
    return(max(x,na.rm=TRUE))
}

min.but.one <- function(x,n=1){
    x[is.na(x)] <- 0
    while (n > 0) {
        x <- x[-which.min(x)]
        n <- n - 1
    }
    return(min(x))
}

col.exc <- -(1:2)
combined.fpkm.wide[,median_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],
                             1,median)]
combined.fpkm.wide[,max_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,max.na.zero)]
combined.fpkm.wide[,all_but_four_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,function(x){min.but.one(x,4)})]
combined.fpkm.wide[,all_but_three_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,function(x){min.but.one(x,3)})]
combined.fpkm.wide[,all_but_two_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,function(x){min.but.one(x,2)})]
combined.fpkm.wide[,all_but_one_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,min.but.one)]
combined.fpkm.wide[,all_expression :=
                       apply(combined.fpkm.wide[,col.exc,with=FALSE],1,min.na.zero)]
combined.fpkm.wide[,all_but_metatherian_expression :=
                        apply(combined.fpkm.wide[,!grepl("^(monodelphis domestica|tree.*|.*_expression)$",
                                                         colnames(combined.fpkm.wide)),
                                                 with=FALSE],1,min.na.zero)]
setkey(combined.fpkm.wide,"tree_genes")

## all genes which are not housekeeping genes which have minimum 
core.placenta.tree.transcriptome.long <- 
    data.table(melt(combined.fpkm.wide[all_but_metatherian_expression >= 10,
                                       ][order(-median_expression)
                                         ],
                    id.vars=c("tree","tree_genes"),
                    variable.name="Species",
                    value.name="FPKM"
                    ))
core.placenta.tree.transcriptome.long[,housekeeping:=FALSE]

core.placenta.tree.transcriptome.long[tree %in% 
                                      genes.to.tree[housekeeping.genes.superset[,tracking_id],
                                                    tree],
                                      housekeeping:=TRUE]
core.placenta.tree.transcriptome.long[,tree:=
                                           reorder(tree,
                                                   combined.fpkm.wide[core.placenta.tree.transcriptome.long[,"tree_genes",with=FALSE],
                                                                      median_expression])]
core.placenta.tree.transcriptome.long[,tree_genes:=factor(tree_genes)]
core.placenta.tree.transcriptome.long[,tree_genes:=
                                           reorder(tree_genes,
                                                   combined.fpkm.wide[core.placenta.tree.transcriptome.long[,"tree_genes",with=FALSE],
                                                                            median_expression])]


core.placenta.tree.transcriptome.genes.shape <-
    combined.fpkm.wide[,list(tree,
                             tree_genes,
                             median_expression,
                             all_but_metatherian_expression,
                             all_expression,
                             all_but_one_expression,
                             all_but_two_expression,
                             all_but_three_expression,
                             all_but_four_expression,
                             max_expression)]

all.trees <- combined.fpkm.wide[!is.na(tree),unique(tree)]
placenta.transcriptome.trees <- 
    core.placenta.tree.transcriptome.long[!is.na(tree),unique(tree)]

pts <- data.table(all.trees=all.trees)
pts <- pts[,placenta.ts:=all.trees %in% placenta.transcriptome.trees]
pts <- pts[,placenta.ts.permissive:=all.trees %in%
                core.placenta.tree.transcriptome.genes.shape[(all_but_four_expression > 10  |
                                                              all_but_metatherian_expression > 10 ) &
                                                             ! is.na(tree),
                                                             unique(tree)]
           ]
pts[,placenta.ts.all:=all.trees %in%
         core.placenta.tree.transcriptome.genes.shape[all_but_metatherian_expression > 10 &
                                                      ! is.na(tree),
                                                      unique(tree)]]
pts[,placenta.ts.one:=all.trees %in%
         core.placenta.tree.transcriptome.genes.shape[all_but_metatherian_expression > 1 &
                                                 ! is.na(tree),
                                                 unique(tree)]]
expressed.housekeeping.genes <-
    core.placenta.tree.transcriptome.genes.shape[max_expression > 10 &
                                            tree %in%
                                            housekeeping.genes.superset[,gene_short_name]
                                            ][,unique(tree)]
pts <- pts[,housekeeping:=all.trees %in% housekeeping.genes.superset[,gene_short_name]]
pts.tree <- pts[,expressed.housekeeping:=all.trees %in% expressed.housekeeping.genes]

save(core.placenta.tree.transcriptome.long,
     core.placenta.tree.transcriptome.genes.shape,
     pts.tree,
     file=args[length(args)])


