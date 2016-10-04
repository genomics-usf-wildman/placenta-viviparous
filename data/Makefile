#!/usr/bin/make -f

SSCONVERT ?=ssconvert
SHELL=/bin/bash
R=R
ROPTS=-q --no-save --no-restore-data
PERL?=perl

space :=
space +=

FPKM_FILES=Supplemental_Table_2_Homo_NEW_UPDATED.txt $(patsubst %.xlsx,%.txt,$(filter-out $(wildcard *Homo*), $(wildcard *_fpkm.xlsx)))


all: combined_fpkm gene_trees housekeeping_genes_superset placenta_core_transcriptome \
	placenta_go_results_bp placenta_go_results_mf hugo_nomenclature gene_annotation \
	mammal_lifehistories.txt placenta_lineage_specific combined_fpkm_per_sample \
	homology_table

mammal_lifehistories.txt:
	wget -O $@ "http://www.esapubs.org/archive/ecol/E084/093/Mammal_lifehistories_v2.txt"

hugo_nomenclature.txt:
	wget -O $@ "http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_ensembl_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"

hugo_nomenclature: hugo_nomenclature.R hugo_nomenclature.txt
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

gene_annotation: ensembl_gene_annotation.R
	$(MODULE) load R/3.2.0; \
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

go_gene_association.txt: go_gene_association.pl go-basic.obo gene_association.goa_human
	$(PERL) $^ |grep -v '^!' > $@

go_gene_association: go_gene_association.R go_gene_association.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

go-basic.obo:
	wget -O $@ "http://purl.obolibrary.org/obo/go/go-basic.obo"

gene_association.goa_human: gene_association.goa_human.gz
	gzip -dc $^ > $@

gene_association.goa_human.gz: 
	wget -O $@ "http://geneontology.org/gene-associations/gene_association.goa_human.gz"

HOM_MouseHumanSequence.rpt: 
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt"

HOM_AllOrganism.rpt:
	wget -O $@ "ftp://ftp.informatics.jax.org/pub/reports/HOM_AllOrganism.rpt"

oma-ensembl.txt.gz:
	wget -O $@ "http://omabrowser.org/All/oma-ensembl.txt.gz"

oma-groups.txt.gz:
	wget -O $@ "http://omabrowser.org/All/oma-groups.txt.gz"

oma-groups.txt: oma-groups.txt.gz
	gzip -dc $^ > $@

oma-ensembl.txt: oma-ensembl.txt.gz
	gzip -dc $^ > $@

oma-groups_long.txt: oma-group-to-long.pl oma-groups.txt oma-ensembl.txt
	$(PERL) $^ > $@

combined_fpkm: combine_fpkm.R HOM_AllOrganism.rpt oma-groups_long.txt fastq_fpkm/combined_fpkms ensembl_orthologs diamond_against_human/gene_to_human_gene.txt gtf_types/gene_id_to_type.txt
	$(R) $(ROPTS) -f $< --args $@ $(wordlist 2,$(words $^),$^)

combined_fpkm_per_sample: combine_fpkm.R HOM_AllOrganism.rpt oma-groups_long.txt fastq_fpkm/combined_fpkms ensembl_orthologs diamond_against_human/gene_to_human_gene.txt gtf_types/gene_id_to_type.txt
	$(R) $(ROPTS) -f $< --args $@ $(wordlist 2,$(words $^),$^)

homology_table: homology_table.R combined_fpkm
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

ensembl_orthologs: ensembl_orthologs.R
	$(R) $(ROPTS) -f $< --args $@

ensembl_orthologs_long: ensembl_orthologs_long.R
	$(R) $(ROPTS) -f $< --args $@

Compara.80.protein.nhx.emf.gz:
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/current_emf/ensembl-compara/homologies/$@" $@

gene_trees: parse_tree.R Compara.80.protein.nhx.emf.gz
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

%.txt: %.xlsx
	$(SSCONVERT) --export-type Gnumeric_stf:stf_assistant -O 'separator="	" quote="" eol="unix"' $< $@

### placenta lineage specific
placenta_lineage_specific: placenta_lineage_specific.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_lineage_specific.txt: ../common_r_code/array_to_text.R placenta_lineage_specific
	$(R) $(ROPTS) -f $< --args placenta_lineage_specific placenta.lineage.specific.analysis $@

### placenta lineage specific
human_mouse_lineage_specific: human_mouse_lineage_specific.R combined_fpkm_per_sample housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

human_mouse_lineage_specific.txt: ../common_r_code/array_to_text.R human_mouse_lineage_specific
	$(R) $(ROPTS) -f $< --args human_mouse_lineage_specific human.mouse.lineage.specific.analysis $@


### placenta specific transcripts and go analysis of them
placenta_core_transcriptome: placenta_core_transcriptome.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_core_transcriptome.txt: ../common_r_code/array_to_text.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args placenta_core_transcriptome core.placenta.transcriptome.long $@

placenta_core_transcriptome_complete: placenta_core_transcriptome_complete.R combined_fpkm \
	housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_core_transcriptome_complete.txt: ../common_r_code/array_to_text.R \
	placenta_core_transcriptome_complete
	$(R) $(ROPTS) -f $< --args placenta_core_transcriptome_complete \
		placenta.core.transcriptome.complete $@

placenta_go_results_mf: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list MF $@

placenta_go_results_bp: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list BP $@

placenta_go_results_permissive_mf: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list.permissive MF $@

placenta_go_results_permissive_bp: go_analysis.R placenta_core_transcriptome
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list.permissive BP $@


placenta_netpath_results: gsea_analysis.R placenta_core_transcriptome \
	Human_NetPath_September_01_2016_Entrezgene.gmt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) placenta.transcriptome.list $(wordlist 3,3,$^) $@

msigdb.v5.1.entrez.gmt:
	curl -c cookiefile -s "http://software.broadinstitute.org/gsea/j_spring_security_check" -d "j_username=foo@example.com" -d "j_password=password"
	curl -b cookiefile "http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.1/msigdb.v5.1.entrez.gmt" -o $@
	rm -f cookefile

kegg.gmt: msigdb.v5.1.entrez.gmt
	grep '^KEGG_' $< > $@

reactome.gmt: msigdb.v5.1.entrez.gmt
	grep '^REACTOME_' $< > $@

biocarta.gmt: msigdb.v5.1.entrez.gmt
	grep '^BIOCARTA_' $< > $@

GSEA_SETS=kegg reactome biocarta netpath
GSEA_SET_ANALYSES=$(foreach gmt,$(GSEA_SETS),placenta_gsea_$(gmt))
$(GSEA_SET_ANALYSES): placenta_gsea_%: %.gmt gsea_analysis.R placenta_core_transcriptome 
	$(R) $(ROPTS) -f $(wordlist 2,2,$^) --args $(wordlist 3,3,$^) placenta.transcriptome.list $< $@


placenta_classification: placenta_classification.R placental_classification.txt
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

placenta_classification_p: placenta_classification_specific.R	\
	combined_fpkm housekeeping_genes_superset					\
	placenta_classification
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

PLACENTA_CLASSIFICATIONS=shape intimacy interdigitation
MF_GO_PLAC_CLASS=$(foreach class, $(PLACENTA_CLASSIFICATIONS),	\
placenta_classification_go_mf_$(class))
$(MF_GO_PLAC_CLASS): placenta_classification_go_mf_%: go_analysis.R	placenta_classification_p
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) 'placenta.classification.significance[["$*"]]' MF $@

.SECONDEXPANSION:

PATHWAY_PLAC_CLASS=$(foreach class, $(PLACENTA_CLASSIFICATIONS),	\
$(foreach gmt,$(GSEA_SETS),placenta_gsea_$(class)_$(gmt)))
pathway_analyses: $(PATHWAY_PLAC_CLASS) $(GSEA_SET_ANALYSES)
$(PATHWAY_PLAC_CLASS): placenta_gsea_%: gsea_analysis.R \
	$$(lastword $$(subst _,$$(space),$$*)).gmt \
	placenta_classification_p
	$(R) $(ROPTS) -f $< --args  $(wordlist 3,3,$^) \
		'placenta.classification.significance[["$(firstword $(subst _,$(space),$*))"]]' \
		$(lastword $(subst _,$(space),$*)).gmt $@


BP_GO_PLAC_CLASS=$(foreach class, $(PLACENTA_CLASSIFICATIONS),	\
placenta_classification_go_bp_$(class))
$(BP_GO_PLAC_CLASS): placenta_classification_go_bp_%: go_analysis.R	placenta_classification_p
	$(R) $(ROPTS) -f $< --args $(wordlist 2,2,$^) 'placenta.classification.significance[["$*"]]' BP $@

ortholog_expression_matrix.txt: ortholog_expression_matrix.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

ortholog_expression_nexus.nex: ortholog_expression_nexus.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@

ortholog_expression_nexus_continuous.nex: ortholog_expression_nexus_continuous.R combined_fpkm housekeeping_genes_superset
	$(R) $(ROPTS) -f $< --args $(wordlist 2,$(words $^),$^) $@
