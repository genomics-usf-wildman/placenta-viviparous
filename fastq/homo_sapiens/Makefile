#!/usr/bin/make -f


bowtie_indexes: hg19.fa
	module load bowtie2/2.1.0; \
	bowtie2-build $< homo_sapiens_hg19_bt2

# Bos_taurus.UMD3.1.dna.toplevel.fa: Bos_taurus.UMD3.1.dna.toplevel.fa.gz
# 	gzip -dc $< > $@
# 
# Bos_taurus.UMD3.1.80.gtf: Bos_taurus.UMD3.1.80.gtf.gz
# 	gzip -dc $< > $@
# 
# Bos_taurus.UMD3.1.dna.toplevel.fa.gz:
#     rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna.toplevel.fa.gz" $@
# 
# Bos_taurus.UMD3.1.80.gtf.gz: 
#	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/gtf/bos_taurus/Bos_taurus.UMD3.1.80.gtf.gz" $@
