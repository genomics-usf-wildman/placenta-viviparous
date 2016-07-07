#!/usr/bin/make -f

SPECIES=spalax_carmeli
ALIGNMENT_SPECIES=nannospalax_galili
GTF=../nannospalax_galili/nannospalax_galili.gtf
FASTA=../nannospalax_galili/nannospalax_galili.toplevel.fa
BOWTIE_INDEX_DIR=./

COMMON_MAKEFILE=../common_makefile

# this is not an ensembl assembly
NOTENSEMBL=1

NREADS=2

# we'll need to use trinity to do the alignment
USE_TRINITY=1
PROTEIN_ALIGNMENT_SPECIES=mus_musculus
PROTEIN_ALIGNMENT_FASTA=Mus_musculus.GRCm38.pep.all.fa

include ../srr_makefile
