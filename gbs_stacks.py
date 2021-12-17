#!/usr/local/env python3
import glob
import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from gbs_methods_stacks import Methods


__author__ = 'duceppemo'
__version__ = 'v0.1'


class GBS(object):
    def __init__(self, args):
        # Command line arguments

        # Analysis type
        self.referenced = args.referenced
        self.de_novo = args.de_novo

        # I/O
        self.input = args.input  # fastq folder
        if args.paired_end:
            self.read_type = 'pe'
        elif args.single_end:
            self.read_type = 'se'
        elif (args.paired_end and args.single_end) \
                or (not args.paired_end and not args.single_end):
            raise Exception('Please choose between "-se" and "-pe"\n')
        self.ion = args.ion_torrent
        self.ref = args.reference  # reference genome
        self.out_folder = args.output  # output folder

        # Metadata
        self.map = args.population_map
        self.barcodes = args.barcodes

        # Performance
        self.cpu = args.threads
        self.parallel = args.parallel

        # Enzymes used
        self.enz1 = args.enzyme1
        self.enz2 = args.enzyme2

        # SNP filtering
        self.min_maf = args.min_maf
        self.min_depth = args.min_depth
        self.max_missing = 1.00 - args.max_missing

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')
        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)

        # Check if folders are not empty
        result = Methods.check_folder_empty(self.input)
        if result == 0:
            raise Exception('Input folder does not contain files with accepted file extensions: {}'.format(
                Methods.accepted_extensions))

        # Get input files and place info in dictionary
        self.sample_dict['raw'] = Methods.get_files(self.input)

        # Check il all samples in population map file are presents as fastq and vice versa
        Methods.check_map(self.sample_dict['raw'], self.map)

        # Check restriction enzyme(s)
        Methods.check_enzymes([self.enz1, self.enz2])

        # Check mode (referenced vs de novo
        Methods.check_mode(self.ref, self.referenced, self.de_novo)

        ############################################################

        # Step completion report files
        done_trimming = self.out_folder + '/done_trimming'
        done_mapping = self.out_folder + '/done_mapping'
        done_ustacks = self.out_folder + '/done_ustacks'
        done_cstacks = self.out_folder + '/done_cstacks'
        done_sstacks = self.out_folder + '/done_sstacks'
        done_tsv2bam = self.out_folder + '/done_tsv2bam'
        done_gstacks = self.out_folder + '/done_gstacks'
        done_populations = self.out_folder + '/done_populations'

        # Output folders to create
        # cleaned = self.out_folder + '/0_cleaned'
        read_length = self.out_folder + '/length_distribution/'
        trimmed = self.out_folder + '/1_trimmed/'
        mapped = self.out_folder + '/2_mapped/'
        stacks = self.out_folder + '/2_stacks/'
        gstacks = self.out_folder + '/3_gstacks/'
        populations = self.out_folder + '/4_populations/'

        # Create ouput folder
        Methods.make_folder(self.out_folder)

        # Trim reads
        if self.ion:  # if Ion Torrent reads
            """
            The raw IonTorent fastq files are typically already demultiplexed and the barcodes and restriction sites 
            trimmed. Those raw fastq files cannot be processed with "process_ragtags". All is done here is quick
            cleaning pass to remove any sequencing adapters missed by the IonServer and to remove low quality
            and short reads.
            """
            print('Processing IonTorrent reads...')

            if not os.path.exists(done_trimming):
                if self.referenced:
                    # Trim reads with bbduk
                    Methods.parallel_trim_reads(Methods.trim_iontorrent, self.sample_dict['raw'], trimmed,
                                                self.cpu, self.parallel)
                else:  # de novo
                    # Find the best read length for trimming
                    df = Methods.parallel_read_length_dist(self.sample_dict['raw'], read_length, self.cpu)
                    # Create plot
                    import plotly.express as px
                    fig = px.line(df, x=df.index, y='Count')
                    fig.write_html(self.out_folder + '/' + 'first_figure.html', auto_open=False)
                    # Detect peak
                    trim_size = Methods.find_peak(df, self.out_folder)

                    # Trim all reads to specific length
                    Methods.parallel_trim_reads(Methods.trim_iontorrent_size, self.sample_dict['raw'], trimmed,
                                                self.cpu, self.parallel, size=trim_size)
                Methods.flag_done(done_trimming)
            else:
                print('Skipping trimming. Already done.')
        else:  # if Illumina reads
            """
            The raw Illumina reads wont typically be demultiplexed by sample and can be cleaned and demultiplexed
            with "process_ragtags" included with Stacks. 
            """
            print('Processing Illumina reads...')
            if self.read_type == 'se':  # if single-end
                # Trim and demultiplex with stacks
                if not os.path.exists(done_trimming):
                    # It's better no to trim the reads because the de novo process needs reads all the same length
                    # Methods.parallel_trim_reads(Methods.trim_illumina_se, self.sample_dict,
                    #                             cleaned, self.cpu, self.parallel)
                    Methods.process_radtags_se(self.input, self.barcodes, trimmed, self.enz1, renz_2=self.enz2)
                    Methods.flag_done(done_trimming)
                else:
                    print('Skipping trimming. Already done.')
            else:  # elif self.read_type == 'pe':  # if paired-end
                # Trim and demultiplex
                if not os.path.exists(done_trimming):
                    # It's better no to trim the reads because the de novo process needs reads all the same length
                    # Methods.parallel_trim_reads(Methods.trim_illumina_pe, self.sample_dict,
                    #                             cleaned, self.cpu, self.parallel)
                    Methods.process_radtags_pe(self.input, self.barcodes, trimmed, self.enz1, renz_2=self.enz2)
                    Methods.flag_done(done_trimming)
                else:
                    print('Skipping trimming. Already done.')

        # Update sample_dict after trimming
        self.sample_dict['trimmed'] = Methods.get_files(trimmed)

        # Map or stack reads
        if self.referenced:  # Map reads
            if self.read_type == 'se':  # if single-end
                if not os.path.exists(done_mapping):
                    Methods.parallel_map_bowtie2_se(mapped, self.ref, self.sample_dict['trimmed'],
                                                    self.cpu, self.parallel)
                    Methods.flag_done(done_mapping)
                else:
                    print('Skipping mapping. Already done.')
            else:  # elif self.read_type == 'pe':  # if paired-end
                if not os.path.exists(done_mapping):
                    Methods.parallel_map_bowtie2_pe(mapped, self.ref, self.sample_dict['trimmed'],
                                                    self.cpu, self.parallel)
                    Methods.flag_done(done_mapping)
                else:
                    print('Skipping mapping. Already done.')
        else:  # elif self.de_novo  # stack reads
            # ustacks
            if not os.path.exists(done_ustacks):
                Methods.parallel_ustacks(self.sample_dict['trimmed'], stacks, self.cpu, self.parallel)
                Methods.flag_done(done_ustacks)
            else:
                print('Skipping ustacks. Already done.')

            # cstacks
            if not os.path.exists(done_cstacks):
                Methods.parallel_cstacks(self.sample_dict['trimmed'], stacks, self.cpu, self.parallel)
                # Methods.cstacks(stacks, self.map, self.cpu)
                Methods.flag_done(done_cstacks)
            else:
                print('Skipping cstacks. Already done.')

            # sstacks
            if not os.path.exists(done_sstacks):
                Methods.parallel_sstacks(self.sample_dict['trimmed'], stacks, self.cpu, self.parallel)
                # Methods.sstacks(stacks, self.map, self.cpu)
                Methods.flag_done(done_sstacks)
            else:
                print('Skipping sstacks. Already done.')

            # tsv2bam
            if not os.path.exists(done_tsv2bam):
                if self.read_type == 'se':  # if single-end
                    Methods.tsv2bam_se(stacks, self.map, self.cpu)
                else:
                    Methods.tsv2bam_pe(stacks, self.map, self.cpu)
                Methods.flag_done(done_tsv2bam)
            else:
                print('Skipping mapping. Already done.')

        # gstacks (call variants)
        if not os.path.exists(done_gstacks):
            print('Calling variants...')
            if self.referenced:
                Methods.call_snps_gstacks(mapped, self.map, gstacks, self.cpu)
            else:
                # Rename bam files
                bam_list = glob.glob(stacks + '*.bam')
                for bam in bam_list:
                    os.rename(bam, bam.replace('.matches', ''))
                Methods.call_snps_gstacks(stacks, self.map, gstacks, self.cpu)
            Methods.flag_done(done_gstacks)
        else:
            print('Skipping variant calling (gstacks). Already done.')

        # Make stats and create SNP VCF file (populations)
        if not os.path.exists(done_populations):
            print('Computing stats...')
            Methods.make_pop_stats(gstacks, populations, self.map, self.cpu)
            Methods.flag_done(done_populations)
        else:
            print('Skipping population statistics (populations). Already done.')

        # Filter
        print('Filtering variants...')
        Methods.vcftools_filter_depth(populations + 'populations.snps.vcf',
                                      populations + 'populations.1.depth_filtered.vcf', self.min_depth)
        Methods.vcftools_filter_maf(populations + 'populations.1.depth_filtered.vcf',
                                    populations + 'populations.2.maf_filtered.vcf', self.min_maf)
        Methods.vcftools_filter_missing(populations + 'populations.2.maf_filtered.vcf',
                                        populations + 'populations.3.missing_filtered.vcf', self.max_missing)
        Methods.ld_filering(populations + 'populations.3.missing_filtered.vcf',
                            populations + 'populations.ld_stats.vcf')
        Methods.vcftools_stats(populations + 'populations.3.missing_filtered.vcf',
                               populations + 'populations.pop_stats.vcf')

        # Filter VCF to only keep homozygous loci
        Methods.filter_out_heterozygous(populations + 'populations.3.missing_filtered.vcf',
                                        populations + 'populations.3.missing_filtered.homo.vcf')

        # Convert VCF to fasta
        tree = populations + '/tree/'
        Methods.make_folder(tree)
        print('Converting VCF to fasta...')
        Methods.vcf2fasta(populations + 'populations.3.missing_filtered.homo.vcf',
                          tree + 'populations.3.missing_filtered.homo.fasta')
        Methods.vcf2fasta(populations + 'populations.3.missing_filtered.vcf',
                          tree + 'populations.3.missing_filtered.fasta')

        # Make tree
        print('Making tree...')
        Methods.make_tree_raxml(tree + 'populations.3.missing_filtered.homo.fasta', tree, self.cpu)
        Methods.make_tree_raxml(tree + 'populations.3.missing_filtered.fasta', tree, self.cpu)


if __name__ == "__main__":
    max_cpu = cpu_count()
    # max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Referenced or de novo GBS pipeline using Stacks')
    parser.add_argument('--referenced',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Dectect SNPs from mapping (bwa->gstacks->populations')
    parser.add_argument('--de-novo',
                        required=False,
                        action='store_true',  # implies default=False
                        help='Dectect SNPs de novo (ustacks->cstacks->gstacks->populations')
    parser.add_argument('-r', '--reference', metavar='/reference_genome.fasta',
                        required=False,
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
    parser.add_argument('--parallel', metavar=2,
                        required=False,
                        type=int, default=2,
                        help='Number of samples to process in parallel.')
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
