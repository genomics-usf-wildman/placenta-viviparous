library("data.table")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
write.table(human.mouse.lineage.specific.analysis[analysis.name=="hum.all" & fdr <= 0.05 & abs(fc) > 2,human_name],
            quote=FALSE,
            sep="",
            col.names=FALSE,
            row.names=FALSE,
            file=args[length(args)])
