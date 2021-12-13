#!/usr/local/env python3
import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from gbs_methods_stacks import Methods
from collections import defaultdict


__author__ = 'duceppemo'
__version__ = 'v0.1'


class GBS(object):
    def __init__(self, args):
        # Command line arguments
        self.input = args.input
        self.out_folder = args.output
        self.cpu = args.threads
        self.ref = args.reference
        if args.paired_end:
            self.read_type = 'pe'
        elif args.single_end:
            self.read_type = 'se'
        elif (args.paired_end and args.single_end)\
                or (not args.paired_end and not args.single_end):
            raise Exception('Please choose between "-se" and "-pe"\n')
        self.ion = args.ion_torrent
        self.map = args.population_map
        self.barcodes = args.barcodes
        self.enz1 = args.enzyme1
        self.enz2 = args.enzyme2
        self.min_maf = args.min_maf
        self.min_depth = args.min_depth
        self.max_missing = 100 - args.max_missing

        # Data
        self.sample_dict = defaultdict(list)

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')
        # Check if number of CPU and memory requested are valid
        self.cpu = Methods.check_cpus(self.cpu)

        # Check if folders are not empty
        result = Methods.check_folder_empty(self.input)
        if result == 0:
            raise Exception('Input folder does not contain files with accepted file extensions: {}'.format(
                Methods.accepted_extensions))

        # Get input files and place info in dictionary
        Methods.get_files(self.input, self.sample_dict)

        # Check il all samples in population map file are presents as fastq and vice versa
        Methods.check_map(self.sample_dict, self.map)

        # Check restriction enzyme(s)
        Methods.check_enzymes([self.enz1, self.enz2])

        # Process and map reads
        trimmed = self.out_folder + '/trimmed/'
        done_trimming = trimmed + '/done_trimming'
        mapped = self.out_folder + '/mapped/'
        done_mapping = mapped + '/done_mapping'
        if self.ion:
            print('Processing IonTorrent reads...')
            # Trim reads with bbduk
            if not os.path.exists(done_trimming):
                Methods.parallel_trim_reads(Methods.trim_iontorrent, self.sample_dict, trimmed, self.cpu)
                Methods.flag_done(done_trimming)
            else:
                print('Skipping trimming. Already done.')

            # Update sample_dict
            self.sample_dict = defaultdict(list)
            Methods.get_files(trimmed, self.sample_dict)

            # Map reads
            if not os.path.exists(done_mapping):
                Methods.parallel_map_bowtie2_se(mapped, self.ref, self.sample_dict, self.cpu)
                Methods.flag_done(done_mapping)
            else:
                print('Skipping mapping. Already done.')
        else:
            print('Processing Illumina reads...')
            if self.read_type == 'se':
                # Trim and demultiplex
                if not os.path.exists(done_trimming):
                    Methods.process_radtags_se(self.input, self.barcodes, trimmed, self.enz1, renz_2=self.enz2)
                    Methods.flag_done(done_trimming)
                else:
                    print('Skipping trimming. Already done.')

                # Update sample dict
                self.sample_dict = defaultdict(list)
                Methods.get_files(trimmed, self.sample_dict)

                # Map
                if not os.path.exists(done_mapping):
                    Methods.parallel_map_bowtie2_se(mapped, self.ref, self.sample_dict, self.cpu)
                    Methods.flag_done(done_mapping)
                else:
                    print('Skipping mapping. Already done.')
            else:  # elif self.read_type == 'pe':
                # Trim and demultiplex
                if not os.path.exists(done_trimming):
                    Methods.process_radtags_pe(self.input, self.barcodes, trimmed, self.enz1, renz_2=self.enz2)
                    Methods.flag_done(done_trimming)
                else:
                    print('Skipping trimming. Already done.')

                # Update sample dict
                self.sample_dict = defaultdict(list)
                Methods.get_files(trimmed, self.sample_dict)

                # Map
                if not os.path.exists(done_mapping):
                    Methods.parallel_map_bowtie2_pe(mapped, self.ref, self.sample_dict, self.cpu)
                    Methods.flag_done(done_mapping)
                else:
                    print('Skipping mapping. Already done.')

        # Call variants
        done_gstacks = self.out_folder + '/done_gstacks'
        if not os.path.exists(done_gstacks):
            print('Calling variants...')
            Methods.call_snps_gstacks(mapped, self.map, self.out_folder, self.cpu)
            Methods.flag_done(done_gstacks)
        else:
            print('Skipping variant calling (gstacks). Already done.')

        # Make stats and create SNP VCF file
        pop_folder = self.out_folder + '/populations/'
        done_populations = self.out_folder + '/done_populations'
        if not os.path.exists(done_populations):
            print('Computing stats...')
            Methods.make_pop_stats(self.out_folder, pop_folder, self.map, self.cpu)
        else:
            print('Skipping population statistics (populations). Already done.')

        # Filter
        print('Filtering variants...')
        Methods.vcftools_filter_depth(pop_folder + 'populations.snps.vcf',
                                      pop_folder + 'populations.depth_filtered.vcf')
        Methods.vcftools_filter_maf(pop_folder + 'populations.depth_filtered.vcf',
                                    pop_folder + 'populations.maf_filtered.vcf')
        Methods.vcftools_filter_missing(pop_folder + 'populations.maf_filtered.vcf',
                                        pop_folder + 'populations.missing_filtered.vcf')
        Methods.ld_filering(pop_folder + 'populations.missing_filtered.vcf',
                            pop_folder + 'populations.ld_stats.vcf')
        Methods.vcftools_stats(pop_folder + 'populations.missing_filtered.vcf',
                               pop_folder + 'populations.pop_stats.vcf')

        # Filter VCF to only keep homozygous loci
        Methods.filter_out_heterozygous(pop_folder + 'populations.missing_filtered.vcf',
                                        pop_folder + 'populations.missing_filtered.homo.vcf')

        # Convert VCF to fasta
        print('Converting VCF to fasta...')
        Methods.vcf2fasta(pop_folder + 'populations.missing_filtered.homo.vcf',
                          pop_folder + 'populations.missing_filtered.homo.fasta')
        Methods.vcf2fasta(pop_folder + 'populations.missing_filtered.vcf',
                          pop_folder + 'populations.missing_filtered.fasta')

        # Make tree
        print('Making tree...')
        Methods.make_tree_raxml(pop_folder + 'populations.missing_filtered.homo.fasta',
                                pop_folder + '/tree_homo', self.cpu)
        Methods.make_tree_raxml(pop_folder + 'populations.missing_filtered.fasta',
                                pop_folder + '/tree', self.cpu)


if __name__ == "__main__":
    max_cpu = cpu_count()
    # max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Reference-based genotyping-by-sequencing (GBS) pipeline using Stacks')
    parser.add_argument('-r', '--reference', metavar='/reference_genome.fasta',
                        required=True,
                        help='Reference genome for read mapping. Mandatory.')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        help='Folder that contains the fastq files. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        help='Folder to hold the result files. Mandatory.')
    parser.add_argument('-e1', '--enzyme1', metavar='mspI',
                        required=True,
                        help='Restriction enzyme used for the GBS procedure. Mandatory.')
    parser.add_argument('-e2', '--enzyme2', metavar='pstI',
                        required=False,
                        help='Restriction enzyme used for the GBS procedure. Optional.')
    parser.add_argument('-p', '--population-map', metavar='/population_map.tsv',
                        required=True,
                        help='A two-column tab-separated file containing a population map of giving samples.'
                             ' See https://catchenlab.life.illinois.edu/stacks/manual/#popmap for details.')
    parser.add_argument('-b', '--barcodes', metavar='/barcodes.tsv',
                        required=True,
                        help='A two-column tab-separated file containing barcode info of giving samples.'
                             ' See https://catchenlab.life.illinois.edu/stacks/manual/#specbc for details.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({}). Optional'.format(max_cpu))
    parser.add_argument('-pe', '--paired-end',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Input fastq files are paired-end. Must chose "-pe" or "-se".')
    parser.add_argument('-se', '--single-end',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Input fastq files are single-end.')
    parser.add_argument('-ion', '--ion-torrent',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Reads are from IonTorrent (different trim parameters). Default is Illumina.')
    parser.add_argument('--min-maf', default=0.05, metavar='0.05',
                        type=float, required=False,
                        help='Minimum allele frequency.')
    parser.add_argument('--min-depth', default=10, metavar='10',
                        type=int, required=False,
                        help='Minimum average depth of coverage.')
    parser.add_argument('--max-missing', default=10, metavar='0.30',
                        type=float, required=False,
                        help='Maximum percentage of missing values.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    GBS(arguments)
