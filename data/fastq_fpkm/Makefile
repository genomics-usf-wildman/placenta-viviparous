#!/usr/bin/make -f

R=R
ROPTS=-q --no-save --no-restore-data


GENE_FPKM_FILES=$(wildcard *_genes.fpkm_tracking)
ISOFORM_FPKM_FILES=$(wildcard *_isoforms.fpkm_tracking)
STAR_LOG_FILES=$(wildcard *_log_star.txt)

SPECIES_WITH_TRINITY_ALIGNMENTS=nannospalax_galili spalax_carmeli ateles_fusciceps
# 
# BAD_ISOFORMS=$(patsubst %,%_isoforms.fpkm_tracking,$(SPECIES_WITH_TRINITY_ALIGNMENTS))
# BAD_GENES=$(patsubst %,%_genes.fpkm_tracking,$(SPECIES_WITH_TRINITY_ALIGNMENTS))
# TRINITY_RSEM_FILES=$(patsubst %,%_trinity_align_rsem_isoforms.txt,$(SPECIES_WITH_TRINITY_ALIGNMENTS))
TRINITY_DIAMOND_FILES=$(patsubst %,%_trinity_diamond.txt,$(SPECIES_WITH_TRINITY_ALIGNMENTS))


all: combined_fpkms

combined_fpkms: combine_fpkms.R protein_id_to_gene_id.txt \
	$(GENE_FPKM_FILES) $(ISOFORM_FPKM_FILES) \
	$(TRINITY_DIAMOND_FILES) \
	$(STAR_LOG_FILES)
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

protein_id_to_gene_id.txt: protein_to_gene.pl \
	../../fastq/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz \
	../../fastq/mouse/Mus_musculus.GRCm38.80.gtf.gz
	./$^ > $@

link_stuff:
	find ../../ -type f \
		\( -iname '*_isoforms.fpkm_tracking' \
		-o -iname '*_isoforms.fpkm_tracking' \
		-o -iname '*_genes.fpkm_tracking' \
		-o -iname '*_trinity_diamond.txt' \
		-o -iname '*_trinity_align_rsem_isoforms.txt' \
        -o -iname 'Log.final.out' \
		\) \
		-not -path '*local*' \
		-print0|xargs -0 ./link_file;
#	rm -f $(BAD_ISOFORMS) $(BAD_GENES)
