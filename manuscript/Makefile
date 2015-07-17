#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: manuscript.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.tex)
	latexmk -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

manuscript.tex: manuscript.Rnw
