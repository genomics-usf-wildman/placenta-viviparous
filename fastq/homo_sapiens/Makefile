#!/usr/bin/make -f

SPECIES=homo_sapiens
GTF=Homo_sapiens.GRCh38.80.gtf
FASTA=Homo_sapiens.GRCh38.dna.toplevel.fa

FASTQ_FILES=101116_lane_5_1.fastq.gz 101116_lane_5_2.fastq.gz

include ../common_makefile

# GSE66622 is Hughes et al. which is a set of human placenta
# comparisons

GSE66622_series_matrix.txt.gz:
	wget -O $@ "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66622/matrix/$@"

GSE66622_sample_info.txt.gz:
	wget -O $@ \
		"ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66622/suppl/GSE66622%5FSampleInfo%2Etxt%2Egz"

## These are SRR from GSE6662 which are from caesarean sections
SRRS=SRR1850945 \
SRR1850946 \
SRR1850951 \
SRR1850952 \
SRR1850953 \
SRR1850954 \
SRR1850957 \
SRR1850958 \
SRR1850965 \
SRR1850966 \
SRR1850975 \
SRR1850976 \
SRR1850983 \
SRR1850984 \
SRR1850987 \
SRR1850988 \
SRR1850991 \
SRR1850992 \
SRR1850993 \
SRR1850994 \
SRR1850995 \
SRR1850996 \
SRR1851003 \
SRR1851004 \
SRR1851011 \
SRR1851012


make_srr_dirs:
	for a in $(SRRS); do \
		mkdir -p $$a; \
		(cd $$a; \
			ln ../../srr_makefile Makefile -sf; \
			echo "SRR=$$a" > srr_num.mk; \
			echo "GTF=../$(GTF)" >> srr_num.mk; \
			echo "SPECIES=$(SPECIES)" >> srr_num.mk; \
			echo "FASTA=../$(FASTA)" >> srr_num.mk; \
		); \
	done;

get_srr: $(patsubst %,%-get_srr,$(SRRS))

$(patsubst %,%-get_srr,$(SRRS)): %-get_srr: %
	make -C $* get_srr;

.PHONY: get_srr $(patsubst %,%-get_srr,$(SRRS))

