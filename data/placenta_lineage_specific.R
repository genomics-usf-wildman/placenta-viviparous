library("data.table")
library("reshape2")
library("org.Hs.eg.db")

args <- c("combined_fpkm","housekeeping_genes_superset","placenta_core_transcriptome")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])

combined.fpkm.wide <- 
    data.table(dcast(combined.fpkm,
                     human_name~species,
                     fun.aggregate=sum,
                     value.var="mean_fpkm"))

### ignore genes which aren't expressed at least 10 FPKM in some
### tissue
placenta.lineage.specific <-
    combined.fpkm.wide[apply(combined.fpkm.wide[,-1,with=FALSE],1,max) >= 100 & !is.na(human_name),]

species <- c("ateles fusciceps", "homo sapiens", "pan paniscus", "mus musculus", 
             "nannospalax galili", "spalax carmeli", "bos taurus", "ovis aries", 
             "sus scrofa", "equus caballus", "canis familiaris", "loxodonta africana", 
             "dasypus novemcinctus", "monodelphis domestica")

### infraclass
eutharians <- species[1:13]
metatheria <- species[14]

### magnorder
boreoeutharians <- species[1:11]
atlantogenata <- species[12:13]

### superorder
euarchontoglires <- species[1:6]
laurasiatheria <- species[7:11]

### order
primates <- species[1:3]
rodentia <- species[4:6]

new.world.monkey <- species[1]
old.world.monkey <- species[2:3]


calculate.t.test <- function(x,g1,g2){
    p.value <- NA;
    try(p.value <- t.test(as.numeric(x[g1]),
                          as.numeric(x[g2]))$p.value,
        silent=TRUE);
    return(p.value)
}


calculate.fc <- function(x,g1,g2){
    fc <- NA
    try({fc <- mean(as.numeric(x[g1]))/mean(as.numeric(x[g2]))
        if (fc < 1) {fc <- -1/fc}
        if (any(!is.finite(fc))) {
            fc <- max(c(mean(as.numeric(x[g1])),mean(as.numeric(x[g2])),na.rm=TRUE))
        }},
        silent=TRUE)
    return(fc)
}
    

calculate.z.test <- function(x,g1,g2){
    p.value <- NA;
    try(p.value <- pnorm(q=as.numeric(x[g1]),
                         mean=mean(as.numeric(x[g2]),na.rm=TRUE),
                         sd=sd(as.numeric(x[g2]),na.rm=TRUE)),
        silent=TRUE);
    return(p.value)
}

comparisons <- list("bor.atl"=list(group1=boreoeutharians,
                                   group1.name="boreoeutharians",
                                   group2=atlantogenata,
                                   group2.name="atlantogenata"
                                   ),
                    "eua.lau"=list(group1=euarchontoglires,
                                   group1.name="euarchontoglires",
                                   group2=laurasiatheria,
                                   group2.name="laurasiatheria"),
                    "rod.pri"=list(group1=rodentia,
                                   group1.name="rodentia",
                                   group2=primates,
                                   group2.name="primates"),
                    "new.old"=list(group1=new.world.monkey,
                                   group1.name="new.world.monkey",
                                   group2=old.world.monkey,
                                   group2.name="old.world.monkey"),
                    "old.all"=list(group1=old.world.monkey,
                                   group1.name="old.world.monkey",
                                   group2=species[!species %in% old.world.monkey],
                                   group2.name="not old.world.monkey"),
                    "rod.all"=list(group1=rodentia,
                                   group1.name="rodentia",
                                   group2=species[!species %in% rodentia],
                                   group2.name="not rodentia"),
                    "pri.all"=list(group1=primates,
                                   group1.name="primates",
                                   group2=species[!species %in% primates],
                                   group2.name="not primates"),
                    "eua.all"=list(group1=euarchontoglires,
                                   group1.name="euarchontoglires",
                                   group2=species[!species %in% euarchontoglires],
                                   group2.name="not euarchontoglires"),
                    "lau.all"=list(group1=laurasiatheria,
                                   group1.name="laurasiatheria",
                                   group2=species[!species %in% laurasiatheria],
                                   group2.name="not laurasiatheria"),
                    "met.eut"=list(group1=metatheria,
                                   group1.name="metatheria",
                                   group2=eutharians,
                                   group2.name="eutharians")
                    )

placenta.lineage.specific.analysis <- list()
for (c.n in names(comparisons)) {
    results <- data.table(placenta.lineage.specific[,human_name])
    setnames(results,"human_name")
    if (length(comparisons[[c.n]]$group1) > 1 &
        length(comparisons[[c.n]]$group2) > 1) {
        results[,p.value:=apply(placenta.lineage.specific,1,
                                      calculate.t.test,
                                      comparisons[[c.n]]$group1,
                                      comparisons[[c.n]]$group2
                                      )]
    } else {
        results[,p.value:=apply(placenta.lineage.specific,1,
                                      calculate.z.test,
                                      comparisons[[c.n]]$group1,
                                      comparisons[[c.n]]$group2
                                      )]
    }
    results[,fdr:=p.adjust(results[,p.value],method="BH")]
    results[,group1.name:=comparisons[[c.n]]$group1.name]
    results[,group2.name:=comparisons[[c.n]]$group2.name]
    results[,fc:=apply(placenta.lineage.specific,1,
                       calculate.fc,
                       comparisons[[c.n]]$group1,
                       comparisons[[c.n]]$group2
                       )
            ]
    results[,group1.mean:=apply(placenta.lineage.specific[,comparisons[[c.n]]$group1,with=FALSE],
                                1,
                                mean,
                                na.rm=TRUE)]
    results[,group1.sd:=apply(placenta.lineage.specific[,comparisons[[c.n]]$group1,with=FALSE],
                              1,
                              sd,
                              na.rm=TRUE)]
    results[,group2.mean:=apply(placenta.lineage.specific[,comparisons[[c.n]]$group2,with=FALSE],
                                1,
                                mean,
                                na.rm=TRUE)]
    results[,group2.sd:=apply(placenta.lineage.specific[,comparisons[[c.n]]$group2,with=FALSE],
                              1,
                              sd,
                              na.rm=TRUE)]
    results[(group2.mean < 50 & group1.mean < 50),p.value:=NA]
    results[,fdr:=p.adjust(results[,p.value],method="BH")]
    results[,analysis.name:=c.n]
    placenta.lineage.specific.analysis[[c.n]] <- results
}
placenta.lineage.specific.analysis <- rbindlist(placenta.lineage.specific.analysis)

save(placenta.lineage.specific.analysis,
     placenta.lineage.specific,
     file=args[length(args)])
