#!/usr/bin/make -f

# this is the env variable to tell us how many processors on this node
# we get
ifdef PBS_NUM_PPN
CORES=$(PBS_NUM_PPN)
else
CORES=1
endif

genes.fpkm_tracking: homo_sapiens_placenta.sorted.bam
	module load cufflinks/2.2.1; \
	cufflinks -p $(CORES) -G gencode.v19.annotation.gtf $<

homo_sapiens_placenta.sorted.bam: homo_sapiens_placenta.bam
	samtools sort $< -o $@ -T temp

homo_sapiens_placenta.bam: 101116_lane_5_1.fastq 101116_lane_5_2.fastq homo_sapiens_hg19_bt2.1.bt2
	module load bowtie2/2.1.0; \
	module load samtools/1.0; \
	bowtie2 -p $(CORES) -x $(wordlist 3,3,$^) \
		-1 $(wordlist 1,1,$^) -2 $(wordlist 2,2,$^) | \
		samtools view -b -o $@ -;

homo_sapiens_hg19_bt2.1.bt2: hg19.fa
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
