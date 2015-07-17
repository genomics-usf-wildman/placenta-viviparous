to.org <- function(x){
    sub("([0-9.-]+)([eE])\\+?(-?)0*([0-9]+)","$\\1\\\\times 10^{\\3\\4}$",
        x)
}
to.latex <- function(x){
  sub("([0-9]+)([eE])\\+?(-?)0*([0-9]+)","\\1\\\\ensuremath{\\\\times 10^{\\3\\4}}",
      x)
}
cleanup.tolatex <- function(output) {
  gsub("\\\\textrm\\{e\\}(-?)0*(\\d+)","\\\\ensuremath{\\\\times 10^{\\1\\2}}",output);
}

trimws <- function(x,left=TRUE,right=TRUE){
  result <- x
  if(left)
    result <- gsub('^\\s+','',x,perl=TRUE)
  if(right)
    result <- gsub('\\s+$','',x,perl=TRUE)

  return(result)
}


print.summary.glm.xtable <- function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...) {

    cat("\\begin{verbatim}\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("Deviance Residuals: \n")
    if (x$df.residual > 5) {
        x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
        names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
            "Max")
    }
    xx <- zapsmall(x$deviance.resid, digits + 1)
    cat("\\end{verbatim}\n")
    print(table.to.latex(xx,digits=digits))
    cat("\\begin{verbatim}\n")
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    }
    else {
        df <- if ("df" %in% names(x)) 
            x[["df"]]
        else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L])) 
            cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
                sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, 
                colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        cat("\\end{verbatim}\n")
        colnames(coefs) <- gsub("(Pr\\(>\\|z\\|\\))","$\\1$",colnames(coefs),perl=TRUE)
        if (nrow(coefs) > 15) {
          print(table.to.latex(coefs,longtable=TRUE))
        } else {
          print(table.to.latex(coefs))
        }
        cat("\\begin{verbatim}\n")
    }
    cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", 
        format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null", 
            "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", 
            "deviance")]), digits = max(5, digits + 1)), " on", 
            format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 
            1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4, digits + 1)), 
        "\n\n", "Number of Fisher Scoring iterations: ", x$iter, 
        "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, 
                  digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\\end{verbatim}\n")
    cat("\n")
    invisible(x)
}



table.to.latex <- function(table,
                           digits=NULL,
                           format=NULL,
                           scientific=NA,
                           nsmall=0,
                           colspec="r",
                           useBooktabs=TRUE,
                           longtable=FALSE,
                           caption=NULL,
                           label=NULL,
                           rownames=TRUE,
                           centering=FALSE,
                           latex.pos="",
                           toprule=if(useBooktabs) "\\toprule" else "\\hline\\hline",
                           midrule=if(useBooktabs) "\\midrule" else "\\hline",
                           bottomrule=if(useBooktabs) "\\bottomrule" else "\\hline\\hline",
                           
                           ...) {
  ans <- ""
  header <- c(if(rownames){""}else{NULL},colnames(table));
  .format.function <- function(x){
    format(if(is.null(format)||is.null(scientific)) {x} else {
      if(!is.na(suppressWarnings(as.numeric(x)))) as.numeric(x) else x},
           digits=digits,nsmall=nsmall,scientific=scientific,...)
  }
  if(is.null(dim(table))){
    header <- names(table)
    rownames <- FALSE
  }
  res <- rbind(header,
               cbind(if(rownames){rownames(table)}else{NULL},
                     if(is.null(dim(table))){t(sapply(table,.format.function))} else
                     {apply(table,1:2,.format.function)}
                     ##format(table,digits=digits,scientific=scientific)
                     ))
  res[,] <- gsub("_","\\\\_",res[,])
  ignore.rows <- -1
  coefrows <- 2:NROW(res)
  if(rownames) {
      coefrows <- 1:NROW(res)
      ignore.rows <- TRUE
  }
##  res[coefrows,ignore.rows] <- sub("(\\*+)","$^{\\1}$",
##                                   res[coefrows,ignore.rows])
  if (NROW(res) > 1) {
      res[coefrows,ignore.rows] <- sub("([0-9]+)([eE])\\+?(-?)0*([0-9]+)","\\1\\\\ensuremath{\\\\times 10^{\\3\\4}}",
                                       res[coefrows,ignore.rows])
  }
  tabspec <- rep("l",ncol(res))
  if (length(colspec)==length(tabspec) &&
      length(colspec) > 1
      ) {
    tabspec <- colspec
  }

  tabbegin <- paste("\\begin{tabular}{",paste(tabspec,collapse=""),"}",sep="")
  tabend <- "\\end{tabular}\n"
  if (longtable) {
    tabbegin <- paste("\\begin{longtable}{",paste(tabspec,collapse=""),"}",sep="")
    tabend <- "\\end{longtable}\n"
    if (!is.null(label)) {
      tabend <- paste(sep="","\\label{",label,"}\\\\","\n",
                      tabend
                      )
    }
    if(!is.null(caption))
      tabend <- paste(sep="","\\caption{",caption,"}\\\\","\n",
                      tabend
                      )
  } else if (!is.null(label)) {
      tabbegin <- paste0("\\begin{table}",latex.pos,"\n",
                         if (centering) {"\\centering\n"} else {""},
                         tabbegin)
      if (!is.null(caption)) {
          tabend <- paste0(tabend,"\\caption{",caption,"}\n")
      }
      tabend <- paste0(tabend,"\\label{",label,"}\n",
                       "\\end{table}","\n")
  }
  ans <- c(ans,tabbegin)
  header <- NULL
  if(length(toprule))
    header <- c(header,toprule)
  if (!rownames)
    header <- c(header,paste("\\multicolumn{1}{c}{",gsub("\\_"," ",trimws(res[1,1])),"}",sep=""))
  for(j in 2:ncol(res))
    header <- c(header,paste(" & \\multicolumn{1}{c}{",gsub("\\_"," ",trimws(res[1,j])),"}",sep=""))
  header <- c(header,"\\\\")
  if(length(midrule))
    header <- c(header,midrule)
  if (longtable)
    header <- c(header,"\\endhead")
  ans <- c(ans,header)
  if (NROW(res) > 1) {
      for(i in 2:NROW(res)) {
          ans <- c(ans,
                   paste(paste(res[i,],collapse=" & "),"\\\\"))
      }
  }
  if(length(bottomrule))
    ans <- c(ans,bottomrule)
  ans <- c(ans,tabend)
  structure(ans,class="Latex")
  
}
