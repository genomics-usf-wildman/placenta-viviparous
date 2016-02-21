#!/usr/bin/make -f

all: uofi_mark.pdf igb_wordmark.pdf relevant_xkcd.png

igb_wordmarks.zip:
	wget -O $@ "http://www.igb.illinois.edu/sites/default/files/upload/Wordmarks-Eps.zip"

igb_wordmark.pdf: igb_wordmarks.zip
	unzip igb_wordmarks.zip Wordmark-black-tag-outlines.eps
	epstopdf --outfile=$@ Wordmark-black-tag-outlines.eps
	rm -f Wordmark-black-tag-outlines.eps

uofi_mark.zip: 
	wget -O $@ "http://identitystandards.illinois.edu/assets/logos/uclogo_1867_horz_bold.zip"

uofi_mark.pdf: uofi_mark.zip
	unzip uofi_mark.zip uclogo_1867_horz_bold/pc/uclogo_1867_horz_bold.eps
	epstopdf --outfile=$@ uclogo_1867_horz_bold/pc/uclogo_1867_horz_bold.eps
	rm -rf uclogo_1867_horz_bold

relevant_xkcd.png:
	wget -O $@ "http://imgs.xkcd.com/comics/dna.png"
