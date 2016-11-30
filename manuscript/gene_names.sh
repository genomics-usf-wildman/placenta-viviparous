#!/bin/bash

for gene in "$@"; do
    name="$(echo "$gene" | gene_to_aliases|awk -F"$(echo -ne "\t")" '{print $2}')"
    echo "\\DeclareAcronym{$gene}{short=\\textit{$gene},long={$name},first-style=reversed,single={\\textit{$gene} ($name)}}";
done;
