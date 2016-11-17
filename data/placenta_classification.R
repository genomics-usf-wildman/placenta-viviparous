library("data.table")


args <- c("../data/placental_classification.txt","placenta_character_states",
          "placental_classification")
args <- commandArgs(trailingOnly=TRUE)


placenta_classification <- fread(args[1])

placenta_classification[,species:=factor(species)]
placenta_classification[,genus:=factor(gsub("^(.)","\\U\\1",gsub(" .*","",species),perl=TRUE))]
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

placenta_classification[,intimacy:=
                             factor(intimacy,
                                    levels=c("epitheliochorial",
                                             "endotheliochorial",
                                             "hemochorial"
                                             ),
                                    ordered=TRUE
                                    )]

### merge in placenta character states from mess/carter 2006
load(args[2])
setkey(placenta_character_states,"Genus")

placenta_classification[,genus:=factor(gsub("^(.)","\\U\\1",gsub(" .*","",species),perl=TRUE))]
setkey(placenta_classification,"genus")
plaenta_classification <-
    placenta_character_states[placenta_classification]
setkey(placenta_classification,
       "species")

save(placenta_classification,
     file=args[length(args)])
