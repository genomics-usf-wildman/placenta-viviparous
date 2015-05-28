#!/usr/bin/make -f


# this is the env variable to tell us how many processors on this node
# we get
ifdef PBS_NUM_PPN
CORES=$(PBS_NUM_PPN)
else
CORES=1
endif

genes.fpkm_tracking: bos_taurus_placenta.sorted.bam Bos_taurus.UMD3.1.80.gtf
	module load cufflinks/2.2.1; \
	cufflinks -p $(CORES) -G $(wordlist 2,2,$^) $<

local/genes.fpkm_tracking: bos_taurus_placenta_local.sorted.bam Bos_taurus.UMD3.1.80.gtf
	module load cufflinks/2.2.1; \
	cufflinks -p $(CORES) -G $(wordlist 2,2,$^) $<

bos_taurus_placenta.sorted.bam: bos_taurus_placenta.bam
	module load samtools/1.0; \
	samtools sort $< -o $@ -T temp

bos_taurus_placenta_local.sorted.bam: bos_taurus_placenta_local.bam
	module load samtools/1.0; \
	samtools sort $< -o $@ -T temp

bos_taurus_placenta.bam: 091125_lane_6_1.fastq 091125_lane_6_2.fastq \
	bos_taurus_umd3_bt2.1.bt2
	module load bowtie2/2.1.0; \
	module load samtools/1.0; \
	bowtie2 -p $(CORES) -x $(patsubst %.1.bt2,%,$(wordlist 3,3,$^)) \
		-1 $(wordlist 1,1,$^) -2 $(wordlist 2,2,$^) | \
		samtools view -b -o $@ -;

bos_taurus_placenta_local.bam: 091125_lane_6_1.fastq 091125_lane_6_2.fastq \
	bos_taurus_umd3_bt2.1.bt2
	module load bowtie2/2.1.0; \
	module load samtools/1.0; \
	bowtie2 -p $(CORES) -N 1 -L 15 -k 10 -x $(patsubst %.1.bt2,%,$(wordlist 3,3,$^)) \
		-1 $(wordlist 1,1,$^) -2 $(wordlist 2,2,$^) | \
		samtools view -b -o $@ -;


bos_taurus_umd3_bt2.1.bt2: Bos_taurus.UMD3.1.dna.toplevel.fa
	module load bowtie2/2.1.0; \
	bowtie2-build $< bos_taurus_umd3_bt2

Bos_taurus.UMD3.1.dna.toplevel.fa: Bos_taurus.UMD3.1.dna.toplevel.fa.gz
	gzip -dc $< > $@

Bos_taurus.UMD3.1.80.gtf: Bos_taurus.UMD3.1.80.gtf.gz
	gzip -dc $< > $@

Bos_taurus.UMD3.1.dna.toplevel.fa.gz:
    rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna.toplevel.fa.gz" $@

Bos_taurus.UMD3.1.80.gtf.gz: 
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-80/gtf/bos_taurus/Bos_taurus.UMD3.1.80.gtf.gz" $@
