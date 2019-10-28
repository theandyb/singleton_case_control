# Exploring the contribution of local sequence context on germline mutation rates using matched case-control sampling

## Introduction

Taking inspiration from Zhu, Neeman, Yap and Huttley (2017) \[[article](https://www.ncbi.nlm.nih.gov/pubmed/27974498)\], this analysis will look at exploring the impact of local sequence context on the germline mutation rate using singleton mutation observations matched with nearby control observations.

## Required Software

It is recommended that this analysis is run in a virtualenv. In particular, I use Anaconda's `conda` to manage my environments and software, with `bioconda` and `conda-forge` channels in addition to the `defaults` channel.

* python 3.6+
* pyfaidx

## Obtaining the Reference Data

### Human Reference Genome (v37)

```bash
mkdir "reference_data/human_g1k_v37"
curl -s "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz" | gunzip -c > "reference_data/human_g1k_v37/human_g1k_v37.fasta"
```

## Directory Structure

* designs: documents detailing code design
* src: source files for analyses
* reference_data: outside data used in analyses
    * human_g1k_v37: v37 of the human reference genome
    
## Running the Analyses

1. `sbatch src/batch_sample.sh` to sample control distribution for singleton observations

        Note: change `src/sampling.py` to sample more/less bases surrounding sampled sites
        
2. `sbatch src/batch_aggregate.sh` to generate motif counts from control observations
3. `sbatch src/batch_singleton_count.sh` to generate motif counts from singletons
