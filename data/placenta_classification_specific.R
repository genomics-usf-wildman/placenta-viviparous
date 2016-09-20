library("data.table")
library("reshape2")
library("org.Hs.eg.db")

args <- c("combined_fpkm","housekeeping_genes_superset",
          "placental_classification.txt","placenta_classification_specific")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
load(args[3])

combined.fpkm.wide <- 
    data.table(dcast(combined.fpkm,
                     human_name~species,
                     fun.aggregate=sum,
                     value.var="mean_fpkm"))

### ignore genes which aren't expressed at least 10 FPKM in some
### tissue and have FPKM >= 1 in more than 1 species
placenta.class.specific <-
    combined.fpkm.wide[apply(combined.fpkm.wide[,-1,with=FALSE],1,max) >= 10 &
                       apply(combined.fpkm.wide[,-1,with=FALSE],1,function(x){sum(x >= 1)}) > 1 &
                       !is.na(human_name),]

### ignore genes which are housekeeping genes
placenta.class.specific <-
    placenta.class.specific[!(human_name %in% housekeeping.genes.superset[,gene_short_name]),]

placenta.class.specific <-
    data.table(melt(placenta.class.specific,id.vars="human_name",
                    variable.name="species",
                    value.name="mean_fpkm"))
setkey(placenta.class.specific,
       "species")
setkey(placenta_classification,"species")
placenta_classification <-
    placenta.class.specific[placenta_classification]      

calculate.aov <- function(gene,factor="barrier") {
    summary(aov(glm(as.formula(paste0("mean_fpkm~",factor,"+0")),
                    data=placenta_classification[human_name==gene,]))
            )[[1]][factor,"Pr(>F)"]
}

placenta.classification.p <-
    data.table(data.frame(shape=sapply(placenta_classification[,unique(human_name)],
                                       calculate.aov,factor="shape"),
                          barrier=sapply(placenta_classification[,unique(human_name)],
                                         calculate.aov,factor="barrier"),
                          interdigitation=sapply(placenta_classification[,unique(human_name)],
                                                 calculate.aov,factor="interdigitation")),
               keep.rownames=TRUE
               )

placenta.classification.p[,shape.fdr:=p.adjust(shape,method="BH")]
placenta.classification.p[,barrier.fdr:=p.adjust(barrier,method="BH")]
placenta.classification.p[,interdigitation.fdr:=p.adjust(interdigitation,method="BH")]
placenta.classification.p[!duplicated(rn),]
placenta.classification.p[,egid:=mget(rn,org.Hs.egSYMBOL2EG,ifnotfound=NA)]
placenta.classification.significance <- list()
for (class in c("shape","barrier","interdigitation")) {
    placenta.classification.significance[[class]] <-
        ifelse(placenta.classification.p[!is.na(egid),paste0(class,".fdr"),with=FALSE] <= 0.05,
               1,0)
    names(placenta.classification.significance[[class]]) <-
        placenta.classification.p[!is.na(egid),egid]
}

save(placenta.classification.p,
     placenta.classification.significance,
     file=args[length(args)])
