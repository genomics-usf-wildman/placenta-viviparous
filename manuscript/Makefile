#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: manuscript.pdf supplemental_material.pdf

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.tex)
	latexmk -f -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

manuscript.tex: manuscript.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp

supplemental_material.tex: supplemental_material.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp
