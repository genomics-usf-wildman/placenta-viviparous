library("data.table")
library("GSEABase")

args <- c("placenta_core_transcriptome","placenta.transcriptome.list",
          "Human_NetPath_September_01_2016_Entrezgene.gmt",
          "placenta_gsea_results_netpath")
args <- commandArgs(trailingOnly=TRUE)

load(args[1])

data <- as.factor(eval(parse(text=args[2])))

gene.lists <- getGmt(args[3],geneIdType=EntrezIdentifier())

fisher_lists <- function(significant,universe,gmt,pathway) {
    in_universe <- length(intersect(universe,unlist(geneIds(gmt))))
    in_pathway <- length(intersect(universe,unlist(geneIds(pathway))))
    sig_in_universe <- length(intersect(significant,unlist(geneIds(gmt))))
    sig_in_pathway <- length(intersect(significant,unlist(geneIds(pathway))))
    fisher.test(matrix(c(sig_in_pathway,in_pathway,
                         sig_in_universe-sig_in_pathway,in_universe-in_pathway),
                       nrow=2,
                       dimnames=list(Significant=c("In","Out"),
                                     Universe=c("In","out")))
                )$p.value
}

p_values <- sapply(names(gene.lists),
                   function(x){fisher_lists(significant=names(data)[data=="1"],
                                            universe=names(data),
                                            gmt=gene.lists,
                                            gene.lists[[x]])})
res.table <- data.table(data.frame(p_values=p_values,fdr=p.adjust(p_values,method="BH")),keep.rownames=TRUE)
setorder(res.table,"p_values")
setnames(res.table,"rn","Pathway")
eval(parse(text=paste0(args[length(args)]," <- res.table")))
eval(parse(text=paste0("save(",args[length(args)],",file=args[length(args)])")))
