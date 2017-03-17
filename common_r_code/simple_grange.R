library("GenomicRanges")
simple.grange <- function(chr,start,end) {
    return(GRanges(seqnames=Rle(paste0("chr",gsub("^chr","",as.character(chr))),
                       rep(1,length(chr))),
                   ranges=IRanges(start=start,
                       end=end)))
}
