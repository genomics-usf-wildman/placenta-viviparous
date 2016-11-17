library("data.table")

args <- c("mess_carter_character_states.txt","placenta_character_states")
args <- commandArgs(trailingOnly=TRUE)

placenta_character_states <- fread(args[1],sep="|")
placenta_character_states[,V1:=NULL]
placenta_character_states[,V12:=NULL]
setnames(placenta_character_states,
         c("Genus",
           "Gross form",
           "Overall structure",
           "Umbilical vessels",
           "Interhemal barrier",
           "Type of trophoblast in barrier",
           "Trophoblastic fenestrae",
           "Hemophagous regions",
           "Areolae",
           "Differentiation of endometrial stroma"))

recode_factor <- function(orig,new.levels) {
    out.fact <- factor(rep_len(NA,length(orig)),
                       levels=unlist(new.levels))
    for (i in names(new.levels)) {
        out.fact[orig==i] <- new.levels[[i]]
    }
    return(out.fact)
}

recode.factors <-
    list("Gross form"=list("1"="diffuse","2"="zonary","3"="cotyledonary",
                           "4"="discoid","5"="double discoid"),
         "Overall structure"=list("1"="labyrinthine","2"="trabecular to villous"),
         "Umbilical vessels"=list("1"="one artery and one vein",
                                  "2"="two arteries and one vein",
                                  "3"="two arteries and two veins"),
         "Interhemal barrier"=list("1"="epitheliochorial","2"="endotheliochorial",
                                   "3"="hemochorial"),
         "Type of trophoblast in barrier"=list("1"="cytotrophoblast",
                                               "2"="syncytiotrophoblast","3"="both"),
         "Trophoblastic fenestrae"=list("1"="absent","2"="present"),
         "Hemophagous regions"=list("1"="present","2"="absent"),
         "Areolae"=list("1"="present","2"="absent"),
         "Differentiation of endometrial stroma"=list("1"="decidual cells atypical or absent",
                                                      "2"="decidual cells present")
         )

for (a in names(recode.factors)) {
    set(placenta_character_states,
        j=which(colnames(placenta_character_states)==a),
        value=recode_factor(placenta_character_states[,a,with=FALSE],recode.factors[[a]]))
}

save(placenta_character_states,
     file=args[length(args)])
