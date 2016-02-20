#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: placenta_viviparous_presentation.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.tex)
	latexmk -f -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

placenta_viviparous_presentation.tex: placenta_viviparous_presentation.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp

version_it: placenta_viviparous_presentation.pdf
	for a in $^; do \
		cp $$a placenta_viviparous_presentation_$${a%%.pdf}_$$(date +%Y%m%d)_$$(git rev-parse --short HEAD).pdf; \
	done;


.PHONY:
