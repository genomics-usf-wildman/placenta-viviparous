#!/usr/bin/make -f

SPECIES=mus_musculus
GTF=Mus_musculus.GRCm38.80.gtf
FASTA=Mus_musculus.GRCm38.dna.toplevel.fa

FASTQ_FILES=101116_lane_7_1.fastq.gz  101116_lane_7_2.fastq.gz

include ../common_makefile


SRRS=SRR243594 \
SRR243595 \
SRR243596 \
SRR243597 \
SRR243598 \
SRR243599 \
SRR243600 \
SRR243601 \
SRR243602 \
SRR243586 \
SRR243587 \
SRR243588 \
SRR243589 \
SRR243590 \
SRR243591 \
SRR243592 \
SRR243593

make_srr_dirs:
	for a in $(SRRS); do \
		mkdir -p $$a; \
		(cd $$a; \
			ln ../../srr_makefile Makefile -sf; \
			echo "SRR=$$a" > srr_num.mk; \
			echo "STAR_OPTIONS=--outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5" >> srr_num.mk; \
		); \
	done;

get_srr: $(patsubst %,%-get_srr,$(SRRS))

$(patsubst %,%-get_srr,$(SRRS)): %-get_srr: %
	make -C $* get_srr;

.PHONY: get_srr $(patsubst %,%-get_srr,$(SRRS))
