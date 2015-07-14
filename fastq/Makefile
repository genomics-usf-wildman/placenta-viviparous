#!/usr/bin/make -f

SUBDIRS=$(dir $(wildcard */Makefile))

remote_files: $(patsubst %/,%-remote_files,$(SUBDIRS))

$(patsubst %/,%-remote_files,$(SUBDIRS)): %-remote_files: %
	make -C $* remote_files

## this file contains a variable that describes how to connect to the
## cluster host where I ran the analyses and will retreive them onto
## the local host for analysis
-include cluster_host.mk


get_results:
	rsync -avmP \
		--include '**_genes.fpkm_tracking' \
		--include '**_isoforms.fpkm_tracking' \
		--include '**trinity_diamond.txt' \
		--include '**_trinity_align_rsem_isoforms.txt' \
        --include '**/' \
	    --exclude '**' \
		$(CLUSTER_HOST)/fastq/ .;

.PHONY: $(patsubst %,%-remote_files,$(SUBDIRS)) remote_files get_results
