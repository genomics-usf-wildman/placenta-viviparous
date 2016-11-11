library("data.table")
library("reshape2")
library("org.Hs.eg.db")
library("parallel")
library("MASS")

args <- c("combined_fpkm","housekeeping_genes_superset",
          "placenta_classification","placenta_classification_specific")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])
load(args[3])

num.cores <- 6

name.to.ensembl <-
    combined.fpkm[species=="homo sapiens",
                  ][!duplicated(gene_short_name),
                    ][,list(gene_short_name,gene_id)]
name.to.ensembl[,egid:=mget(gene_short_name,org.Hs.egSYMBOL2EG,ifnotfound=NA)]

setkey(name.to.ensembl,
       "gene_short_name")

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

set.seed(30621)
### replace 0s with very small non-zero to avoid trouble with
### underestimating the variance
placenta_classification[mean_fpkm==0,
                        mean_fpkm:=rnorm(sum(mean_fpkm==0,na.rm=TRUE))/100]

calculate.aov <- function(gene,factor="intimacy") {
    temp.data <- placenta_classification[human_name==gene,]
    temp.glm <- glm(as.formula(paste0("mean_fpkm~",factor,"+0")),
                    data=placenta_classification[human_name==gene,])
    temp.glm.coef <- summary(temp.glm)$coefficients
    temp.aov.summary <-
        summary(aov(temp.glm))[[1]]
    rbindlist(list(data.table("gene"=gene,
                              "egid"=name.to.ensembl[gene,egid],
                              "gene_id"=name.to.ensembl[gene,gene_id],
                              "factor"=factor,
                              "type"="aov",
                              "level"=factor,
                              "coefficient"=temp.aov.summary[factor,"Mean Sq"],
                              "statistic"=temp.aov.summary[factor,"F value"],
                              "p"=temp.aov.summary[factor,"Pr(>F)"]),
                   data.table("gene"=gene,
                              "egid"=name.to.ensembl[gene,egid],
                              "gene_id"=name.to.ensembl[gene,gene_id],
                              "factor"=factor,
                              "type"="glm",
                              "level"=gsub(factor,"",names(temp.glm.coef[,"t value"])),
                              "coefficient"=temp.glm.coef[,"Estimate"],
                              "statistic"=temp.glm.coef[,"t value"],
                              "p"=temp.glm.coef[,"Pr(>|t|)"])))
}

calculate.polr <- function(gene,factor="intimacy") {
    temp.coef <- data.frame("Value"=NA,"t value"=NA,check.names=FALSE)
    rownames(temp.coef) <- "mean_fpkm"
    try({
        temp.coef <-
            coef(summary(polr(as.formula(paste0(factor,"~mean_fpkm")),
                              placenta_classification[human_name==gene],
                              Hess=TRUE)))},
        silent=TRUE)
    data.table("gene"=gene,
               "egid"=name.to.ensembl[gene,egid],
               "gene_id"=name.to.ensembl[gene,gene_id],
               "factor"=factor,
               "type"="polr",
               "level"=factor,
               "coefficient"=temp.coef["mean_fpkm","Value"],
               "statistic"=temp.coef["mean_fpkm","t value"],
               ## approximate p using the Z distribution
               "p"=pnorm(abs(as.numeric(temp.coef["mean_fpkm","t value"])),
                         lower.tail=FALSE)*2)
}
    

placenta.classification.p <-
    rbindlist(lapply(c("shape","intimacy","interdigitation"),
                     function(x){
                         rbindlist(mclapply(placenta_classification[,unique(human_name)],
                                            calculate.aov,factor=x,mc.cores=num.cores))})
              )
placenta.classification.p <-
    rbindlist(c(list(orig=placenta.classification.p),
                lapply(c("intimacy","interdigitation"),
                       function(x){
                           rbindlist(mclapply(placenta_classification[,unique(human_name)],
                                              calculate.polr,factor=x,mc.cores=num.cores))})
                ))

## throw out analyses with infinite statistic
placenta.classification.p <-
    placenta.classification.p[is.finite(statistic),]
    


## adjust for fdr by type and class
placenta.classification.p[,fdr:=p.adjust(p,method="BH"),by=list(factor,type)]
placenta.classification.p[,fdr.overall:=p.adjust(p,method="BH")]

placenta.classification.significance <- list()
for (class in c("shape","intimacy","interdigitation")) {
    placenta.classification.significance[[class]] <-
        ifelse(placenta.classification.p[!is.na(egid) & type=="aov" & factor==class,
                                         fdr] <= 0.05,
               1,0)
    names(placenta.classification.significance[[class]]) <-
        placenta.classification.p[!is.na(egid) & type=="aov" & factor==class,egid]
}

save(placenta.classification.p,
     placenta.classification.significance,
     file=args[length(args)])
