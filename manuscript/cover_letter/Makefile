#!/usr/bin/make -f

EMACSARGS=--batch  -l ~/.emacs 
EMACS=emacs

ORG_FILES?=$(wildcard *.org)
ORG_FILES=cover_letter.org

PDF_FROM_ORG_FILES=$(patsubst %.org,%.pdf,$(ORG_FILES))

DOC_FROM_ORG_FILES=$(patsubst %.org,%.doc,$(ORG_FILES))

all: pdf

pdf: $(PDF_FROM_ORG_FILES)

doc: $(DOC_FROM_ORG_FILES)

uofi_mark.zip: 
	wget -O $@ "http://identitystandards.illinois.edu/assets/logos/uclogo_1867_horz_bold.zip"

uofi_mark.pdf: uofi_mark.zip
	unzip uofi_mark.zip uclogo_1867_horz_bold/pc/uclogo_1867_horz_bold.eps
	epstopdf --outfile=$@ uclogo_1867_horz_bold/pc/uclogo_1867_horz_bold.eps
	rm -rf uclogo_1867_horz_bold

$(PDF_FROM_ORG_FILES): %.pdf: %.org uofi_mark.pdf
	$(EMACS) $(EMACSARGS) --visit $< --funcall org-latex-export-to-pdf

$(DOC_FROM_ORG_FILES): %.doc: %.org uofi_mark.pdf
	$(EMACS) $(EMACSARGS) --visit $< --eval '(setq org-odt-preferred-output-format "doc")' --funcall org-odt-export-to-odt

.PHONY: doc pdf
