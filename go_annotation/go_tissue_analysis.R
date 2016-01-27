library("data.table")
library(doMC)
ncore = multicore:::detectCores()
registerDoMC(cores = ncore)
library("topGO")
library("org.Hs.eg.db")

args <- c("combined_fpkm","housekeeping_genes_superset","go_bos_taurus_mf")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
load(args[2])



do.go.analysis <- function(species.name,go.ontology,top.100=TRUE,exclude.housekeeping=FALSE) {

    all.genes <- combined.fpkm[species==species.name & !is.na(human_name),]
    all.genes[,egid:=as.vector(unlist(sapply(mget(all.genes[,human_name],
                   org.Hs.egSYMBOL2EG,ifnotfound=NA),function(x){x[1]})))]
    all.genes <- all.genes[!is.na(egid),]

    analysis.type <- ""
    if (exclude.housekeeping) {
        all.genes <- all.genes[!(human_name %in% housekeeping.genes.superset[,gene_short_name]),]
        analysis.type <- " No Housekeeping"
    }
    
    all.genes[,rank:=rank(-mean_fpkm)]

    if (top.100) {
        analysis.type <- paste0("Top 100",analysis.type)
        gene.list <- all.genes[,factor(as.integer(rank <=100))]
    } else {
        analysis.type <- paste0("Top 1%",analysis.type)
        gene.list <- all.genes[,factor(as.integer(rank <= (sum(mean_fpkm>1)*0.01)))]
    }
    names(gene.list) <- all.genes[,egid]
    
    godata <- new("topGOdata",
                  ontology=go.ontology,
                  allGenes=gene.list,
                  annot=annFUN.org,
                  mapping="org.Hs.eg.db")
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    result.fisher <- getSigGroups(godata,test.stat)
    res.table <- 
        data.table(GenTable(godata,
                            classic=result.fisher,
                            orderBy="classic",
                            topNodes=length(result.fisher@score)))
    res.table[,classic:=as.numeric(gsub("< ","",classic))]
    setnames(res.table,"classic","p value")
    res.table[,FDR:=p.adjust(method="BH",`p value`)]
    res.table[,Ontology:=go.ontology]
    res.table[,species:=species.name]
    res.table[,type:=analysis.type]
    return(res.table)
}
 
results <- 
    c(lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="MF",
             top.100=TRUE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="MF",
             top.100=TRUE,
             exclude.housekeeping=TRUE
             ),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="MF",
             top.100=FALSE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="MF",
             top.100=FALSE,
             exclude.housekeeping=TRUE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="BP",
             top.100=TRUE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="BP",
             top.100=TRUE,
             exclude.housekeeping=TRUE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="BP",
             top.100=FALSE),
      lapply(combined.fpkm[,unique(species)],
             do.go.analysis,
             go.ontology="BP",
             top.100=FALSE,
             exclude.housekeeping=TRUE))

print(results)

go.results <- rbindlist(results)

go.results <- go.results[`p value` <= 0.05,]

save(file=args[length(args)],
     go.results)
               
