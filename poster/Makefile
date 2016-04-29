#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: references.bib placenta_viviparous_poster.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.cls)
	latexmk -f -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

references_clean.bib: placenta_viviparous_poster.tex references.bib
	biber placenta_viviparous_poster --output-format bibtex --output-file $@

references.bib:
	if ! [ -e references.bib ]; then \
		ln references_clean.bib $@; \
	else \
		touch $@; \
	fi;

placenta_viviparous_poster.pdf: placenta_viviparous_poster.tex acronyms.tex

placenta_viviparous_poster.tex: placenta_viviparous_poster.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp

version_it: placenta_viviparous_poster.pdf
	for a in $^; do \
		cp $$a placenta_viviparous_$${a%%.pdf}_$$(date +%Y%m%d)_$$(git rev-parse --short HEAD).pdf; \
	done;


.PHONY:
