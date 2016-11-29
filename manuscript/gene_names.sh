#!/bin/bash

for gene in "$@"; do
    name="$(echo "$gene" | gene_to_aliases|awk -F"$(echo -ne "\t")" '{print $2}')"
    echo "\\DeclareAcronym{$gene}{short=\\textit{$gene},long={$name}}";
done;
