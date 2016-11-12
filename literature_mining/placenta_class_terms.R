library("data.table")

args <- commandArgs(trailingOnly=TRUE)

load(args[1])
write.table(c(placenta_classification[,levels(shape)],
              placenta_classification[,levels(interdigitation)],
              placenta_classification[,levels(intimacy)]),
            quote=FALSE,
            sep="",
            col.names=FALSE,
            row.names=FALSE,
            file=args[length(args)])
