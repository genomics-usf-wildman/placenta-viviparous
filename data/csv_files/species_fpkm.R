library("data.table")
args <- c("../combined_fpkm","homo_sapiens","homo_sapiens_mean_fpkm.txt")
args <- commandArgs(trailingOnly=TRUE)
this.species <- gsub("_"," ",args[2])

load(args[1])
combined.fpkm[,`:=`(FPKM=NULL,
                    FPKM_conf_lo=NULL,
                    FPKM_conf_hi=NULL,
                    FPKM_status=NULL
                    )]
write.table(file=args[length(args)],
            combined.fpkm[species==this.species,][order(-mean_fpkm),],
            sep="\t")

