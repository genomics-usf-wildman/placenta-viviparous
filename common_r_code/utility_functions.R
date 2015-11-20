data.matrix.dt <- function(x){
    temp <- data.matrix(as.data.frame(x)[,-1])
    rownames(temp) <- as.data.frame(x)[,1]
    temp
}

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                             {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

capfirst <- function(s, strict = FALSE) {
    paste(toupper(substring(s, 1, 1)),
          {s <- substring(s, 2); if(strict) tolower(s) else s},
          sep = "")
}

species.names <- 
    list(
        "ateles fusciceps"=
            list(sort.order=1,
                 common.name="Spider monkey",
                 scientific.name="Ateles fuciceps"),
        "homo sapiens"=
            list(sort.order=2,
                 common.name="Human",
                 scientific.name="Homo sapiens"),
        "pan paniscus"=
            list(sort.order=3,
                 common.name="Bonobo",
                 scientific.name="Pan paniscus"),
        "mus musculus"=
            list(sort.order=4,
                 common.name="mouse",
                 scientific.name="Mus musculus"),
        "nannospalax galili"=
            list(sort.order=5,
                 name="nannospalax galili",
                 scientific.name="Nannospalax galili"),
        "spalax carmeli"=
            list(sort.order=6,
                 name="spalax carmeli",
                 scientific.name="Spalax carmeli"),
        "bos taurus"=
            list(sort.order=7,
                 name="Cow",
                 scientific.name="Bos taurus"),
        "ovis aries"=
            list(sort.order=8,
                 name="Sheep",
                 scientific.name="Ovis aries"),
        "sus scrofa"=
            list(sort.order=9,
                 name="Boar",
                 scientific.name="Sus scrofa"),
        "equus caballus"=
            list(sort.order=10,
                 name="Horse",
                 scientific.name="Equus caballus"),
        "canis familiaris"=
            list(sort.order=11,
                 name="Dog",
                 scientific.name="Canis familiaris"),
        "loxodonta africana"=
            list(sort.order=12,
                 name="Elephant",
                 scientific.name="Loxodonta africana"),
        "dasypus novemcinctus"=
            list(sort.order=13,
                 name="Nine banded armadillo",
                 scientific.name="Dasypus novemcinctus"),
        "monodelphis domestica"=
            list(sort.order=14,
                 name="Opposum",
                 scientific.name="Monodelphis domestica"))


species <- function(s,join=TRUE) {
    temp <- paste0("\\textit{",capfirst(s),"}")
    if (join) {
        return(english.join(temp))
    } else {
        return(temp)
    }
}

english.join <- function(s){
    ret <- c()
    ## everything but the last gets a , appended if there are more than two
    if (length(s) > 2) {
        s[-length(s)] <- 
            paste0(s[-length(s)],",")
    }
    ## the second to last gets an and appended if there are more than one
    if (length(s) > 1) {
        s[length(s)-1] <- 
            paste0(s[length(s)-1]," and")
    }
    return(paste(sep="",collapse=" ",s))
}


source("/home/don/projects/xtable/pkg/R/print.xtable.R")

