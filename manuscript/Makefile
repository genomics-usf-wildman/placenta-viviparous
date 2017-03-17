#!/usr/bin/make -f


R ?= R
ROPTS=-q --no-save --no-restore-data

all: references.bib manuscript.pdf supplemental_material.pdf \
	extended_discussion.pdf extended_supplemental_material.pdf

bmcart.cls:
	wget -O temp.zip "http://media.biomedcentral.com/content/production/bmc_article-tex.zip"
	unzip -j temp.zip '*/bmcart.cls';
	rm -f temp.zip

%.pdf: %.svg
	inkscape -D -A $@ $<

%.png: %.svg
	inkscape -D -e $@ -d 300 $<

## we require double evaluation here because we use load_cache to
## replay the results of analyses
%.tex: %.Rnw
	$(R) --encoding=utf-8 -e "library('knitr'); knit('$<')"
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

TEX_DEPS=common_latex.tex acronyms.tex

manuscript.pdf: manuscript.tex $(TEX_DEPS) supplemental_material.aux

extended_discussion.pdf: extended_discussion.tex $(TEX_DEPS)	\
supplemental_material.aux

supplemental_material.pdf: manuscript.tex $(TEX_DEPS)

extended_supplemental_material.pdf:				\
extended_supplemental_material.tex $(TEX_DEPS)

GO_GSEA_RESULTS= $(wildcard ../data/placenta_go_results_??)			\
$(wildcard ../data/placenta_go_results_permissive_??) $(wildcard	\
../data/placenta_gsea_*) $(wildcard									\
../data/placenta_classification_go_??)

manuscript.tex: manuscript.Rnw ../data/combined_fpkm	\
../data/placenta_core_transcriptome						\
../data/placenta_core_tree_transcriptome				\
../data/placenta_core_tree_transcriptome_permissive		\
../data/human_mouse_lineage_specific					\
../data/placenta_lineage_specific $(GO_GSEA_RESULTS)

supplemental_material.tex: supplemental_material.Rnw	\
../data/combined_fpkm ../data/placenta_lineage_specific	\
../data/human_mouse_lineage_specific					\
../data/placenta_core_tree_transcriptome				\
../data/placenta_core_tree_transcriptome_permissive		\
../data/placenta_core_transcriptome						\
../data/orthologs_by_species_ranking					\
../data/placenta_core_transcriptome $(GO_GSEA_RESULTS)

extended_supplemental_material.tex: extended_supplemental_material.Rnw	\
../data/combined_fpkm ../data/placenta_lineage_specific	\
../data/human_mouse_lineage_specific					\
../data/placenta_core_transcriptome $(GO_GSEA_RESULTS)

version_it: manuscript.pdf supplemental_material.pdf			\
	extended_supplemental_material.pdf extended_discussion.pdf
	for a in $^; do \
		cp $$a placenta_viviparous_$${a%%.pdf}_$$(date +%Y%m%d)_$$(git rev-parse --short HEAD).pdf; \
	done;


.PHONY:
