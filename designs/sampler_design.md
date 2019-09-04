# Design Document: Sampling Control Distribution

## Problem
Given a file listing observations of singletons (tab-delimited) and a reference genome (fasta), match each observed singleton with a control observation from the reference genome.

## Input Structure

### Singletons

|CHR|POS|REF|ALT|AA|AN|Motif|Category|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|20|60155|A|G|.|7522|CTTTAAGAG(CTCTTAAAG)|AT_GC|

### Reference Genome

Fasta file with chromosomes numbered without prefix or suffix (i.e. '22', not 'chr22').

## Procedure

1. Create Fasta object pointing to reference fasta file
2. Read singletons one-by-one
    * Get position and reference allele
    * Get sequence from fasta in 300bp window around position
    * Get indices of location of reference allele in string, ignoring center allele
    * Sample 1 index, return 9-mer sequence centered at site