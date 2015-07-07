#!/usr/bin/make -f

SPECIES=nannospalax_galili
ALIGNMENT_SPECIES=nannospalax_galili
GTF=nannospalax_galili.gtf
FASTA=nannospalax_galili.toplevel.fa

COMMON_MAKEFILE=../common_makefile
# this is not an ensembl assembly
NOTENSEMBL=1

include ../srr_makefile


# get the fasta files; the perl regex is to just keep the ref name
# which is used in the gtf file
nannospalax_galili.toplevel.fa:
	rsync -avP \
	"rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/CHR_Un/1026970_ref_S.galili_v1.0_chrUn.fa.gz" \
	"rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/CHR_MT/1026970_ref_S.galili_v1.0_chrMT.fa.gz" \
	.
	gzip -dc 1026970_ref_S.galili_v1.0_chrMT.fa.gz 1026970_ref_S.galili_v1.0_chrUn.fa.gz | \
		perl -pe 's/^>gi\|\d+\|ref\|([^\|]+)\|.+/>$$1/; $$_;' > \
		$@

# build the gtf file from the gff3 file
nannospalax_galili.gtf: nannospalax_galili.toplevel.fa
	rsync -avP "rsync://ftp.ncbi.nih.gov/genomes/Nannospalax_galili/GFF/ref_S.galili_v1.0_top_level.gff3.gz" .
	gzip -dc ref_S.galili_v1.0_top_level.gff3.gz > ref_S.galili_v1.0_top_level.gff3
	$(MODULE) load cufflinks/2.2.1; \
	gffread -o nannospalax_galili.gtf -T -g nannospalax_galili.toplevel.fa ref_S.galili_v1.0_top_level.gff3
