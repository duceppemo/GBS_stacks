# GBS_stacks
A pipeline to process Genotyping-by-Sequencing (GBS) data from Illumina or IonTorrent (single-end or paired-end) using `stacks` to call the genotypes and the variants.

## Installation
Requirements are install throuh conda:
```
# Create environment (requires conda installed)
conda create -n gbs_stacks -c bioconda stacks raxml bowtie2 fastp bbmap pysam

# Activate environment
conda activate gbs_stacks

# Clone repository (requires git installed)
git clone https://github.com/duceppemo/GBS_stacks

# Test installation
cd GBS_stacks
python gbs_stacks.py -h
```

## Usage
```
usage: gbs_stacks.py [-h] -r /reference_genome.fasta -i /input_folder/ -o /output_folder/ -e1 mspI [-e2 pstI] -p /population_map.tsv -b /barcodes.tsv [-t 48] [-pe] [-se] [-ion]

Reference-based genotyping-by-sequencing (GBS) pipeline using Stacks

optional arguments:
  -h, --help            show this help message and exit
  -r /reference_genome.fasta, --reference /reference_genome.fasta
                        Reference genome for read mapping. Mandatory.
  -i /input_folder/, --input /input_folder/
                        Folder that contains the fastq files. Mandatory.
  -o /output_folder/, --output /output_folder/
                        Folder to hold the result files. Mandatory.
  -e1 mspI, --enzyme1 mspI
                        Restriction enzyme used for the GBS procedure. Mandatory.
  -e2 pstI, --enzyme2 pstI
                        Restriction enzyme used for the GBS procedure. Optional.
  -p /population_map.tsv, --population-map /population_map.tsv
                        A two-column tab-separated file containing a population map of giving samples.
                        See https://catchenlab.life.illinois.edu/stacks/manual/#popmap for details.
  -b /barcodes.tsv, --barcodes /barcodes.tsv
                        A two-column tab-separated file containing barcode info of giving samples.
                        See https://catchenlab.life.illinois.edu/stacks/manual/#specbc for details.
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48). Optional
  -pe, --paired-end     Input fastq files are paired-end. Must chose "-pe" or "-se".
  -se, --single-end     Input fastq files are single-end.
  -ion, --ion-torrent   Reads are from IonTorrent (different trim parameters). Default is Illumina.
```
