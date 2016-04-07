#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: references.bib manuscript.pdf supplemental_material.pdf

bmcart.cls:
	wget -O temp.zip "http://media.biomedcentral.com/content/production/bmc_article-tex.zip"
	unzip -j temp.zip '*/bmcart.cls';
	rm -f temp.zip

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"

%.pdf: %.tex $(wildcard *.bib) $(wildcard *.cls)
	latexmk -f -pdf -pdflatex='xelatex -interaction=nonstopmode %O %S' -bibtex -use-make $<

references_clean.bib: manuscript.tex references.bib
	biber manuscript --output-format bibtex --output-file $@

references.bib:
	if ! [ -e references.bib ]; then \
		ln references_clean.bib $@; \
	else \
		touch $@; \
	fi;

manuscript.pdf: manuscript.tex common_latex.tex acronyms.tex supplemental_material.aux

supplemental_material.pdf: manuscript.tex common_latex.tex acronyms.tex

manuscript.tex: manuscript.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp

supplemental_material.tex: supplemental_material.Rnw ../data/combined_fpkm ../data/placenta_core_transcriptome ../data/placenta_go_results_mf ../data/placenta_go_results_bp

version_it: manuscript.pdf supplemental_material.pdf
	for a in $^; do \
		cp $$a placenta_viviparous_$${a%%.pdf}_$$(date +%Y%m%d)_$$(git rev-parse --short HEAD).pdf; \
	done;


.PHONY:
