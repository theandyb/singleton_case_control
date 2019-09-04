# Exploring the contribution of local sequence context on germline mutation rates using matched case-control sampling

## Introduction

Taking inspiration from Zhu, Neeman, Yap and Huttley (2017) \[[article](https://www.ncbi.nlm.nih.gov/pubmed/27974498)\], this analysis will look at exploring the impact of local sequence context on the germline mutation rate using singleton mutation observations matched with nearby control observations.

## Required Software

It is recommended that this analysis is run in a virtualenv. In particular, I use Anaconda's `conda` to manage my environments and software, with `bioconda` and `conda-forge` channels in addition to the `defaults` channel.

* python 3.6+
* pyfaidx

## Directory Structure

* designs: documents detailing code design
* src: source files for analyses
* reference_data: 