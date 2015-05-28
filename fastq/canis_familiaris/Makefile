#!/usr/bin/make -f

# this is the env variable to tell us how many processors on this node
# we get
ifdef PBS_NUM_PPN
CORES=$(PBS_NUM_PPN)
else
CORES=1
endif

genes.fpkm_tracking: canis_familiaris_placenta.sorted.bam canis_familiaris_broadd2.67.gtf
	module load cufflinks/2.2.1; \
	cufflinks -p $(CORES) -G $(wordlist 2,2,$^) $<

canis_familiaris_placenta.sorted.bam: canis_familiaris_placenta.bam
	module load samtools/1.0; \
	samtools sort $< -o $@ -T temp

canis_familiaris_placenta.bam: 101116_lane_8_1.fastq 101116_lane_8_2.fastq canis_familiaris_broadd2_bt2.1.bt2
	module load bowtie2/2.1.0; \
	module load samtools/1.0; \
	bowtie2 -p $(CORES) -x $(patsubst %.1.bt2,%,$(wordlist 3,3,$^)) \
		-1 $(wordlist 1,1,$^) -2 $(wordlist 2,2,$^) | \
		samtools view -b -o $@ -;

canis_familiaris_broadd2_bt2.1.bt2: canis_familiaris_broadd2.67.dna.toplevel.fa
	module load bowtie2/2.1.0; \
	bowtie2-build $< canis_familiaris_broadd2_bt2

canis_familiaris_broadd2.67.dna.toplevel.fa: canis_familiaris_broadd2.67.dna.toplevel.fa.gz
	gzip -dc $< > $@

canis_familiaris_broadd2.67.gtf: canis_familiaris_broadd2.67.gtf.gz
	gzip -dc $< > $@

canis_familiaris_broadd2.67.dna.toplevel.fa.gz:
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-67/fasta/canis_familiaris/dna/Canis_familiaris.BROADD2.67.dna_rm.toplevel.fa.gz" $@

canis_familiaris_broadd2.67.gtf.gz: 
	rsync -avP "rsync://ftp.ensembl.org/ensembl/pub/release-67/gtf/canis_familiaris/Canis_familiaris.BROADD2.67.gtf.gz" $@
