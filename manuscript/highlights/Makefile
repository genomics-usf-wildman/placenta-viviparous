#!/usr/bin/make -f

EMACSARGS=--batch  -l ~/.emacs 
EMACS=emacs

ORG_FILES?=$(wildcard *.org)
ORG_FILES=highlights.org

PDF_FROM_ORG_FILES=$(patsubst %.org,%.pdf,$(ORG_FILES))

DOC_FROM_ORG_FILES=$(patsubst %.org,%.doc,$(ORG_FILES))

all: pdf

pdf: $(PDF_FROM_ORG_FILES)

doc: $(DOC_FROM_ORG_FILES)

$(PDF_FROM_ORG_FILES): %.pdf: %.org
	$(EMACS) $(EMACSARGS) --visit $< --funcall org-latex-export-to-pdf

$(DOC_FROM_ORG_FILES): %.doc: %.org
	$(EMACS) $(EMACSARGS) --visit $< --eval '(setq org-odt-preferred-output-format "doc")' --funcall org-odt-export-to-odt

.PHONY: doc pdf
