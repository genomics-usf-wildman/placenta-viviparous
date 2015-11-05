#!/usr/bin/make -f


placenta_core_genes.txt: ../data/placenta_core_transcriptome.txt
	awk '{print $$2}' $<|sort|uniq|sed 's/"//g'|grep -v Species > $@

core_vs_terms.txt: pubmed_grid placenta_core_genes.txt placenta_terms.txt
	touch core_vs_terms_pmids.txt
	./$< --terms $(wordlist 2,2,$^) --terms $(wordlist 3,3,$^) \
		--pmid-input core_vs_terms_pmids.txt \
		--pmid-output core_vs_terms_pmids.txt > $@

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
