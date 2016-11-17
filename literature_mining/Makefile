#!/usr/bin/make -f

R=R
ROPTS=-q --no-save --no-restore-data

placenta_class_genes.txt: placenta_class_specific.R ../data/placenta_classification_p
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_class_terms.txt: placenta_class_terms.R		\
	../data/placenta_classification placenta_terms.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) $@
	cat $(wordlist 3,3,$^) >> $@

human_lineage_genes.txt: human_lineage_genes.R ../data/human_mouse_lineage_specific
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

mouse_lineage_genes.txt: mouse_lineage_genes.R ../data/human_mouse_lineage_specific
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_core_genes.txt: ../data/placenta_core_transcriptome.txt
	awk '{print $$2}' $<|sort|uniq|sed 's/"//g'|grep -v Species > $@

class_vs_terms.txt: pubmed_grid placenta_class_genes.txt placenta_class_terms.txt
	touch class_vs_terms_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 3,3,$^) \
		--pmid-input class_vs_terms_pmids.txt \
		--pmid-output class_vs_terms_pmids.txt > $@

core_vs_terms.txt: pubmed_grid placenta_core_genes.txt placenta_terms.txt
	touch core_vs_terms_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 3,3,$^) \
		--pmid-input core_vs_terms_pmids.txt \
		--pmid-output core_vs_terms_pmids.txt > $@

core_vs_core.txt: pubmed_grid placenta_core_genes.txt
	touch core_vs_core_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 2,2,$^) \
		--pmid-input core_vs_core_pmids.txt \
		--pmid-output core_vs_core_pmids.txt > $@

human_lineage_vs_terms.txt: pubmed_grid human_lineage_genes.txt placenta_terms.txt
	touch human_lineage_vs_terms_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 3,3,$^) \
		--pmid-input human_lineage_vs_terms_pmids.txt \
		--pmid-output human_lineage_vs_terms_pmids.txt > $@

mouse_lineage_vs_terms.txt: pubmed_grid mouse_lineage_genes.txt placenta_terms.txt
	touch mouse_lineage_vs_terms_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 3,3,$^) \
		--pmid-input mouse_lineage_vs_terms_pmids.txt \
		--pmid-output mouse_lineage_vs_terms_pmids.txt > $@

core_only.txt: pubmed_grid placenta_core_genes.txt
	touch core_only_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --dimension 1 \
		--pmid-input core_only_pmids.txt \
		--pmid-output core_only_pmids.txt > $@

terms_only.txt: pubmed_grid placenta_terms.txt
	touch terms_only_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --dimension 1 \
		--pmid-input terms_only_pmids.txt \
		--pmid-output terms_only_pmids.txt > $@
