#!/usr/bin/make -f

SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data

GENE_TREES=cga csh1 cyp19a1 

ALIGNMENTS=$(patsubst %,%_gene_tree_genes_promoter_only.aln,$(GENE_TREES))

DOTPLOT_ALIGNMENTS=$(patsubst %,%_gene_tree_promoter_alignment.pdf,$(GENE_TREES))
DOTPLOT_FASTA=$(patsubst %,%_gene_tree_promoter_dotplot.pdf,$(GENE_TREES))

FIGURES=$(DOTPLOT_ALIGNMENTS) $(DOTPLOT_FASTA)

all: $(ALIGNMENTS) $(FIGURES) \
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

$(DOTPLOT_ALIGNMENT): %_gene_tree_promoter_alignment.pdf: \
	%_promoter_alignment.R %_gene_tree_genes_promoter_only.aln ../dotplots.R
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

$(DOTPLOT_FASTA): %_gene_tree_promoter_dotplot.pdf: \
	%_promoter_dotplot.R %_gene_tree_genes_promoter_only.fasta ../dotplots.R
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

cyp19a1_gene_tree_promoter_alignment.pdf: cyp19a1_promoter_alignment.R cyp19a1_gene_tree_genes_promoter_only.aln ../dotplots.R
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
