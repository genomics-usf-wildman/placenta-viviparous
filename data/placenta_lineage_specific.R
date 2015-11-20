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
evolution.of.expression <-
    combined.fpkm.wide[apply(combined.fpkm.wide[,-1,with=FALSE],1,max) >= 100,]

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


calculate.t.test <- function(x,g1,g2){
    p.value <- NA;
    try(p.value <- t.test(as.numeric(x[g1]),
                          as.numeric(x[g2]))$p.value,
        silent=TRUE);
    return(p.value)
}

calculate.z.test <- function(x,g1,g2){
    p.value <- NA;
    try(p.value <- pnorm(q=as.numeric(x[g1]),
                         mean=mean(as.numeric(x[g2]),na.rm=TRUE),
                         sd=sd(as.numeric(x[g2]),na.rm=TRUE)),
        silent=TRUE);
    return(p.value)
}

evolution.of.expression[,bor.atl.p:=apply(evolution.of.expression,1,
                             calculate.t.test,boreoeutharians,atlantogenata)]
evolution.of.expression[,bor.atl.fdr:=p.adjust(bor.atl.p,method="BH")]

evolution.of.expression[,eua.lau.p:=apply(evolution.of.expression,1,
                             calculate.t.test,euarchontoglires,laurasiatheria)]
evolution.of.expression[,eua.lau.fdr:=p.adjust(eua.lau.p,method="BH")]

evolution.of.expression[,rod.pri.p:=apply(evolution.of.expression,1,
                             calculate.t.test,rodentia,primates)]
evolution.of.expression[,rod.pri.fdr:=p.adjust(rod.pri.p,method="BH")]

evolution.of.expression[,met.eut.p:=apply(evolution.of.expression,1,
                             calculate.z.test,metatheria,eutharians)]
evolution.of.expression[,met.eut.fdr:=p.adjust(met.eut.p,method="BH")]

save(evolution.of.expression,
     file=args[length(args)])
