#!/usr/bin/make -f

SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data

GENE_TREES=cga csh1 cyp19a1 

ALIGNMENTS=$(patsubst %,%_gene_tree_genes_promoter_only.aln,$(GENE_TREES))

all: $(ALIGNMENTS) \
	$(patsubst %,%_gene_tree_genes.fasta,$(GENE_TREES))

$(ALIGNMENTS): %.aln: %.fasta
	clustalw $<

csh1_gene_tree_genes.txt: get_genes_in_tree.R ../../data/combined_fpkm ../../data/gene_trees
	$(R) $(ROPTS) -f $< --args CSH1 $@ $(wordlist 2,$(words $^),$^)

%_gene_tree_genes.fasta: %_gene_tree_genes.txt
	../get_gene_sequences.pl $< > $@

%_gene_tree_genes_promoter_only.fasta: %_gene_tree_genes.txt
	../get_gene_sequences.pl --flank 5000 --promoter-only $< > $@

%_gene_tree_genes_only.fasta: %_gene_tree_genes.txt
	../get_gene_sequences.pl --flank 0 $< > $@

cga_gene_tree_genes.txt: get_genes_in_tree.R ../../data/combined_fpkm ../../data/gene_trees
	$(R) $(ROPTS) -f $< --args CGA $@ $(wordlist 2,$(words $^),$^)

cyp19a1_gene_tree_genes.txt: get_genes_in_tree.R ../../data/combined_fpkm ../../data/gene_trees
	$(R) $(ROPTS) -f $< --args CYP19A1 $@ $(wordlist 2,$(words $^),$^)

csh1_gene_tree_promoter_alignment.pdf: csh1_promoter_alignment.R csh1_gene_tree_genes_promoter_only.aln ../dotplots.R
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
