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
conda create -n gbs_stacks -c bioconda -c ploly stacks raxml bowtie2 fastp bbmap pysam git \
    psutil pandas matplotlib scipy numpy plotly ete3

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
usage: python gbs_stacks.py [-h] [--referenced] [--de-novo] [-r /reference_genome.fasta] -i /input_folder/ -o /output_folder/ [-s 64|auto] -e1 mspI [-e2 pstI] -p /population_map.tsv -b /barcodes.tsv [-t 48] [--parallel 2] [-pe] [-se] [-ion]
                            [--min-maf 0.05] [--min-depth 10] [--max-missing 0.30]

Referenced or de novo GBS pipeline using Stacks

optional arguments:
  -h, --help            show this help message and exit
  --referenced          Dectect SNPs from mapping (bwa->gstacks->populations
  --de-novo             Dectect SNPs de novo (ustacks->cstacks->gstacks->populations
  -r /reference_genome.fasta, --reference /reference_genome.fasta
                        Reference genome for read mapping. Mandatory.
  -i /input_folder/, --input /input_folder/
                        Folder that contains the fastq files. Mandatory.
  -o /output_folder/, --output /output_folder/
                        Folder to hold the result files. Mandatory.
  -s 64|auto, --size 64|auto
                        Minimum read size to keep after trimming. An experimetial auto size selection (use "auto" as argument) is also available. It is based on highest peak detection after plotting read length distribution.
                        Experimental. Default is 64. Optional.
  -e1 mspI, --enzyme1 mspI
                        Restriction enzyme used for the GBS procedure. Mandatory.
  -e2 pstI, --enzyme2 pstI
                        Restriction enzyme used for the GBS procedure. Optional.
  -p /population_map.tsv, --population-map /population_map.tsv
                        A two-column tab-separated file containing a population map of giving samples. See https://catchenlab.life.illinois.edu/stacks/manual/#popmap for details. Mandatory.
  -b /barcodes.tsv, --barcodes /barcodes.tsv
                        A two-column tab-separated file containing barcode info of giving samples. See https://catchenlab.life.illinois.edu/stacks/manual/#specbc for details. Mandatory.
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48). Optional
  --parallel 2          Number of samples to process in parallel.
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
    -i ~/input_fastq_folder \
    -se \
    -ion \
    -o ~/output_gbs_folder \
    -e1 pstI \
    -e2 mspI
```
You must change the paths found in the help or in this example to your own paths. Just copy/pasting this command will not work on your system.
## Inputs
Here are example of the required inputs:
* Barcode description (barcode.tsv):
```commandline
CTAAGGTAA	sample1
TAAGGAGAA	sample2
```
* Population map (pop_map_stacks.tsv):
```commandline
sample1	sample1
sample2	sample2
```
* Absolute path of folder location containing the fastq files (folder structure from "tree" command):
```
/input_fastq_folder
├── sample1_R1.fastq.gz
└── sample2_R1.fastq.gz
```
## Outputs
See `stacks` manual for output files description (https://catchenlab.life.illinois.edu/stacks/manual/)

Other output files of interest:
* /output_gbs_folder/populations/populations.snps.vcf -> unfiltered SNPs from gstacks
* /output_gbs_folder/populations/populations.depth_filtered.vcf -> SNPs filtered for minimum average coverage
* /output_gbs_folder/populations/populations.maf_filtered.vcf -> SNPs filted for minimum average coverage + Minimum Allele Frequency
* /output_gbs_folder/populations/populations.missing_filtered.vcf -> SNPs filtered for minimum average coverage + Minimum Allele Frequency + maximum missing values
* /output_gbs_folder/populations/populations.missing_filtered.homo.vcf -> Only sites that display homozygous alleles
* /output_gbs_folder/populations/tree/RAxML_bestTree.populations.missing_filtered.tree -> Maximum likelyhood tree of sample based on SNPs
## SNP filtering
The pipeline has a "resume" function implemented. For example, if the pipeline crashes at some point because an input file is not formatted correctly, the concerned file can be edited and the pipeline re-launched with the same command line and it will resume where it crashed. This is done by writing "report" files once a step is complete.

As a bonus, this implementation allows testing multiple SNP filtering conditions. All that needs to be done is re-run the pipeline by changing the `--min-maf`, `--min-depth` and/or `--max-missing` value(s). Note that the previously generated files will be overwritten doing such. Make sure to move the files to another location before re-launching the pipeline if you want to compare the effect of various filtering settings.
