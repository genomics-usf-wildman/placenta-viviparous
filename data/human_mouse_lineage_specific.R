library("data.table")
library("reshape2")
library("org.Hs.eg.db")
library("lme4")

args <- c("combined_fpkm_per_sample","housekeeping_genes_superset","human_mouse_lineage_specific")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

species <- c("ateles fusciceps", "homo sapiens", "pan paniscus", "mus musculus", 
             "spalax galili", "spalax carmeli", "bos taurus", "ovis aries", 
             "sus scrofa", "equus caballus", "canis familiaris", "loxodonta africana", 
             "dasypus novemcinctus", "monodelphis domestica")


placenta.lineage.specific <-
    data.table(human_name=combined.fpkm[FPKM >= 100 & !is.na(human_name),unique(human_name)])

### ignore genes which are housekeeping genes
placenta.lineage.specific <-
    placenta.lineage.specific[!(human_name %in% housekeeping.genes.superset[,gene_short_name]),]
setkey(combined.fpkm,"human_name")

calculate.mann.whitney.test <- function(x,g1,g2){
    data.subset <- combined.fpkm[x,list(species,FPKM)]
    data.subset[,this.group:=NA]
    data.subset[species %in% g1,this.group:=TRUE]
    data.subset[species %in% g2,this.group:=FALSE]
    data.subset <- data.subset[!is.na(this.group)]
    p.value <- NA;
    try({
        this.sum <-
            summary(lmer(FPKM~1+this.group+(1|species),
                         data.subset))
        p.value <- pt(this.sum$coefficients["this.groupTRUE",
                                            "t value"],
                      df=this.sum$ngrps,lower.tail=FALSE)
    },silent=TRUE);
    return(p.value)
}


calculate.fc <- function(x,g1,g2){
    data.subset <- combined.fpkm[x]
    fc <- NA
    try({
        first <- combined.fpkm[x][species %in% g1,FPKM]
        second <- combined.fpkm[x][species %in% g2,FPKM]
        fc <- mean(as.numeric(first))/mean(as.numeric(second))
        if (fc < 1) {fc <- -1/fc}
        if (any(!is.finite(fc))) {
            fc <- max(c(mean(as.numeric(first)),mean(as.numeric(second)),na.rm=TRUE))
        }},
        silent=TRUE)
    return(fc)
}

comparisons <- list("hum.all"=list(group1="homo sapiens",
                                   group1.name="homo sapiens",
                                   group2=species[!species %in% "homo sapiens"],
                                   group2.name="not humans"
                                   ),
                    "mou.all"=list(group1="mus musculus",
                                   group1.name="mus musculus",
                                   group2=species[!species %in% "mus musculus"],
                                   group2.name="not mice"
                                   )
                    )

human.mouse.linage.specific.analysis <- list()
for (c.n in names(comparisons)) {
    results <- data.table(placenta.lineage.specific[,human_name])
    setnames(results,"human_name")
    results[,p.value:=apply(placenta.lineage.specific,1,
                            calculate.mann.whitney.test,
                            comparisons[[c.n]]$group1,
                            comparisons[[c.n]]$group2
                            )]
    results[,fdr:=p.adjust(results[,p.value],method="BH")]
    results[,group1.name:=comparisons[[c.n]]$group1.name]
    results[,group2.name:=comparisons[[c.n]]$group2.name]
    results[,fc:=apply(placenta.lineage.specific,1,
                       calculate.fc,
                       comparisons[[c.n]]$group1,
                       comparisons[[c.n]]$group2
                       )
            ]
    results[,fdr:=p.adjust(results[,p.value],method="BH")]
    results[,analysis.name:=c.n]
    human.mouse.linage.specific.analysis[[c.n]] <- results
}
human.mouse.linage.specific.analysis <- rbindlist(human.mouse.linage.specific.analysis)

save(human.mouse.linage.specific.analysis,
     placenta.lineage.specific,
     file=args[length(args)])
