# GBS_stacks
A pipeline to process Genotyping-by-Sequencing (GBS) data from Illumina or IonTorrent (single-end or paired-end) using `stacks` to call the genotypes and the variants.
## Installation
To run this pipeline you need a Linux computer. This pipeline has been developed and tested on Ubuntu 20.04 LTS.
### Conda
You need conda installed to run this pipeline. If you already have conda installed, you can skip this step. I recommend using `miniconda`. Instructions to install from the terminal:
```
# Download the latest Linux 64-bit version of miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Allow install file to be executed
chmod +x Miniconda3-latest-Linux-x86_64.sh

# Launch the install script
./Miniconda3-latest-Linux-x86_64.sh
```
Follow the prompt to complete the installation. I recommend accepting (saying "yes") to adding conda auto-activation to your PATH. Once the installation completed, you should close your terminal and open a new one. You should now see `base` at the start of the command line.
### Pipeline
Pipeline requirements are installed through conda:
```
# Create a new conda environment that will hold all the dependencies
conda create -n gbs_stacks -c bioconda stacks raxml bowtie2 fastp bbmap pysam git

# Activate environment
conda activate gbs_stacks

# Install the pipeline (clone repository)
git clone https://github.com/duceppemo/GBS_stacks

# Test installation
cd GBS_stacks
python gbs_stacks.py -h
```
A help message should display without errors if the conda environment was set up properly and the git repository was cloned successfully.
## Usage
```
usage: python gbs_stacks.py [-h] -r /reference_genome.fasta -i /input_folder/ -o /output_folder/ -e1 mspI [-e2 pstI] -p /population_map.tsv -b /barcodes.tsv [-t 48] [-pe] [-se] [-ion] [--min-maf 0.05] [--min-depth 10] [--max-missing 0.30]

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
                        A two-column tab-separated file containing a population map of giving samples. See https://catchenlab.life.illinois.edu/stacks/manual/#popmap for details.
  -b /barcodes.tsv, --barcodes /barcodes.tsv
                        A two-column tab-separated file containing barcode info of giving samples. See https://catchenlab.life.illinois.edu/stacks/manual/#specbc for details.
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48). Optional
  -pe, --paired-end     Input fastq files are paired-end. Must chose "-pe" or "-se".
  -se, --single-end     Input fastq files are single-end.
  -ion, --ion-torrent   Reads are from IonTorrent (different trim parameters). Default is Illumina.
  --min-maf 0.05        Minimum allele frequency.
  --min-depth 10        Minimum average depth of coverage.
  --max-missing 0.30    Maximum percentage of missing values.
```
Here's an example to run the pipeline with several fastq files generated in an IonTorrent sequencer (single-end):
```
# Activate environment, if not already done
conda activate gbs

# Run the pipeline
python gbs_stacks.py -r Secale_cereale_Weining.fasta \
    -b barcodes.tsv \
    -p pop_map_stacks.tsv \
    -i /input_fastq_folder \
    -se \
    -ion \
    -o /output_gbs_folder \
    -e1 pstI \
    -e2 mspI
```
