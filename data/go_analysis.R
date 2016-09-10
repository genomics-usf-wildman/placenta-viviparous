library("data.table")
library("topGO")
library("org.Hs.eg.db")

args <- c("placenta_specific","placenta.transcriptome.list",
          "MF","placenta_go_results_mf")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])

data <- as.factor(eval(parse(text=args[2])))
godata <- new("topGOdata",
              ontology=args[3],
              allGenes=data,
              annot=annFUN.org,
              mapping="org.Hs.eg.db")
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
result.fisher <- getSigGroups(godata,test.stat)
res.table <- 
    data.table(GenTable(godata,
                        classic=result.fisher,
                        orderBy="classic",
                        topNodes=length(result.fisher@score)))
res.table[,classic:=gsub("< ","",classic)]
res.table[,classic.fdr:=p.adjust(method="BH",classic)]
res.table[,Ontology:=args[3]]

eval(parse(text=paste0(args[4]," <- res.table")))

eval(parse(text=paste0("save(",args[4],",file=args[4])")))
