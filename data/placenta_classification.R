library("data.table")


args <- c("../data/placental_classification.txt","placental_classification")
args <- commandArgs(trailingOnly=TRUE)

placenta_classification <- fread(args[1])

placenta_classification[,species:=factor(species)]
placenta_classification[,shape:=factor(shape)]
placenta_classification[,interdigitation:=
                             factor(interdigitation,
                                    levels=c("villous",
                                             "folded",
                                             "trabecular",
                                             "lamellar",
                                             "labyrinthine"
                                             ),
                                    ordered=TRUE
                                    )]

placenta_classification[,barrier:=
                             factor(barrier,
                                    levels=c("epitheliochorial",
                                             "endotheliochorial",
                                             "hemochorial"
                                             ),
                                    ordered=TRUE
                                    )]

save(placenta_classification,
     file=args[length(args)])
