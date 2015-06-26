
dotplot.grid <- function(seq1,seq2,wsize = 1, wstep = 1, nmatch=1) {
    mkwin <- function(t.seq, wsize, wstep) {
        sapply(seq(from = 1, to = length(t.seq) - wsize + 1, by = wstep), 
               function(i) c2s(t.seq[i:(i + wsize - 1)]))
    }
    wseq1 <- mkwin(seq1, wsize, wstep)
    wseq2 <- mkwin(seq2, wsize, wstep)
    "%==%" <- function(x, y) {
        colSums(sapply(x, s2c) ==
                    sapply(y, s2c))
    }
    xy <- outer(wseq1, wseq2, "%==%")
    xy <- xy/nmatch*100
    xy[xy>100] <- 100
    temp <- xy
    temp[] <- heat.colors(101)[round(xy)+1]
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
            dotplot.grid(s2c(as.character(aln[i]))[start:end],s2c(as.character(aln[j]))[start:end],
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
