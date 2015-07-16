
dotplot.grid <- function(seq1,seq2,wsize = 1, wstep = 1, nmatch=wsize) {
    if (missing(nmatch)) {
        nmatch <- wsize
    }
    if (wsize > 1) {
        ## the way this was originally done was totally the most
        ## inefficient way possible
        equality <- outer(seq1,seq2,"==")
        seq1.seq <-
            seq(from = 1, to = length(seq1) - wsize + 1, by = wstep)
        seq2.seq <-
            seq(from = 1, to = length(seq2) - wsize + 1, by = wstep)
        xy <- matrix(NA,nrow=length(seq1.seq),ncol=length(seq2.seq))
        ## now, with a window, average over the whole thing
        for (i in 1:length(seq1.seq)) {
            for (j in 1:length(seq2.seq)) {
                xy[i,j] <- sum(equality[matrix(c(seq1.seq[i]:(seq1.seq[i]+wsize-1),
                                                 seq2.seq[j]:(seq2.seq[j]+wsize-1)
                                                 ),
                                               dimnames=list(NULL,c("row","col")),
                                               ncol=2
                                               )])
            }
        }
    } else {
        xy <- outer(seq1,seq2,"==")
    }
    xy <- xy/nmatch*100-25
    xy[xy>75] <- 75
    xy[xy<0] <- 0
    temp <- xy
    temp[] <- grey.colors(76,start=0,end=1)[76-round(xy)]
    # pushViewport(plotViewport(c(5,4,4,2)))
#     pushViewport(dataViewport(xData = NULL, yData = NULL, xscale = c(1,length(wseq1)),
#                               yscale = c(1,length(wseq2)), extension = 0))
    grid.raster(temp,width=unit(1,"npc"),height=unit(1,"npc"))
    grid.rect()
    # grid.xaxis(name = "xa")
    # grid.yaxis(name = "ya")
    # popViewport()
#    popViewport()

}

dotplot.all <- function(alignment,start,end,seqs=TRUE,...) {
    aln <- alignment$seq[seqs]
    seq.names <- alignment$nam[seqs]
    pushViewport(plotViewport(c(0,7,12,0)))
    pushViewport(viewport(width=min(convertY(unit(1,"npc"),"inches"),convertX(unit(1,"npc"),"inches")),
                          height=min(convertY(unit(1,"npc"),"inches"),convertX(unit(1,"npc"),"inches"))))
    pushViewport(viewport(layout=grid.layout(nrow=length(aln),ncol=length(aln))))
    if (missing(start))
        start <- 1
    if (missing(end)) {
        end <- length(s2c(as.character(aln[1])))
    }
    for (i in 1:length(aln)) {
        for (j in 1:length(aln)) {
            if (i == j)
                next;
            if (j < i)
                next;
            pushViewport(viewport(layout.pos.row=i,layout.pos.col=j))
            dotplot.grid(s2c(as.character(aln[i]))[start:end],
                         s2c(as.character(aln[j]))[start:end],
                         ...
                         )
            popViewport()
        }
    }
    popViewport()
    grid.text(c("",seq.names[-1]),
              x=unit(seq(0.5/length(seq.names),
                  1-0.5/length(seq.names),
                  length.out=length(seq.names)),"npc"),
              y=unit(1,"npc")+unit(0.5,"lines"),
              just="left",
              rot=90)
    grid.text(c(seq.names[-length(seq.names)],""),
              y=unit(rev(seq(0.5/length(seq.names),
                  1-0.5/length(seq.names),
                  length.out=length(seq.names))),"npc"),
              x=unit(seq(1/length(seq.names),
                  1,
                  length.out=length(seq.names)),"npc")-unit(0.5,"lines"),
              just="right")
    popViewport(2)
}
