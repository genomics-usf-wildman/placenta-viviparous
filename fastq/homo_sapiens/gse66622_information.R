library(data.table)

sample.info <- fread("bash -c \"echo -ne '\t'; zcat GSE66622_sample_info.txt.gz\"")
setkey(sample.info,sample)

series.matrix <-
    fread("zgrep -e Sample_title -e geo -e relation GSE66622_series_matrix.txt.gz",
          sep="\t")
series.matrix <- data.frame(t(series.matrix),stringsAsFactors=FALSE)
colnames(series.matrix) <- gsub("!","",series.matrix[1,])
series.matrix <- series.matrix[-1,]
series.matrix$sample <- rownames(series.matrix)
series.matrix <- data.table(series.matrix)
setkey(series.matrix,sample)

series.matrix <-
    series.matrix[sample.info]
series.matrix[,SRX:=gsub(".+(SRX\\d+)","\\1",Sample_relation.1)]
### cat out the subset which are cesarean
cat(paste(collapse="\n",series.matrix[birth=="cesarean",SRX]),"\n")

### to figure out the srr list, output this through srx_to_srr.pl
