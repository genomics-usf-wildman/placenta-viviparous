#!/usr/bin/make -f

SUBDIRS=$(dir $(wildcard */Makefile))

remote_files: $(patsubst %/,%-remote_files,$(SUBDIRS))

$(patsubst %/,%-remote_files,$(SUBDIRS)): %-remote_files: %
	make -C $* remote_files

.PHONY: $(patsubst %,%-remote_files,$(SUBDIRS)) remote_files
