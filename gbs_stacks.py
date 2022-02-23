#!/usr/local/env python3
import glob
import os
from argparse import ArgumentParser
from multiprocessing import cpu_count
from gbs_methods_stacks import Methods
import plotly.express as px


__author__ = 'duceppemo'
__version__ = 'v0.1'


class GBS(object):
    def __init__(self, args):
        # Analysis type
        self.referenced = args.referenced
        self.de_novo = args.de_novo

        # I/O
        self.input = os.path.abspath(args.input)  # fastq folder
        if args.paired_end:
            self.read_type = 'pe'
        elif args.single_end:
            self.read_type = 'se'
        elif (args.paired_end and args.single_end) \
                or (not args.paired_end and not args.single_end):
            raise Exception('Please choose between "-se" and "-pe"\n')
        self.ion = args.ion_torrent
        self.ref = args.reference  # reference genome
        self.out_folder = os.path.abspath(args.output)  # output folder

        # Read length auto removal
        self.size = args.size

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
        self.max_missing = round(1 - args.max_missing, 2)

        # Data
        self.sample_dict = dict()

        # Run
        self.run()

    def run(self):
        print('Checking a few things...')
        # Check if number of CPU and memory requested are valid
        self.cpu, self.parallel = Methods.check_cpus(self.cpu, self.parallel)

        # Check if folders are not empty
        # Convert input to absolute path if relative
        result, self.input = Methods.check_input_folder(self.input)
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
        read_length = self.out_folder + '/0_length_distribution/'
        trimmed = self.out_folder + '/1_trimmed/'
        mapped = self.out_folder + '/2_mapped/'
        stacks = self.out_folder + '/2_stacks/'
        gstacks = self.out_folder + '/3_gstacks/'
        populations = self.out_folder + '/4_populations/'

        # Create ouput folder
        Methods.make_folder(self.out_folder)

        # Auto size selection
        if self.size == 'auto':
            # Find the best read length for trimming
            df = Methods.parallel_read_length_dist(self.sample_dict['raw'], self.cpu)

            # Create plot
            Methods.make_folder(read_length)
            fig1 = px.line(df, x=df.index, y='Count')
            with open(read_length + '/' + 'first_figure.html', 'w') as fh:
                fh.write(fig1.to_html(full_html=False, include_plotlyjs='cdn'))

            # Detect peak
            self.size = Methods.find_peak(df, read_length)

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
                Methods.parallel_trim_reads(Methods.trim_iontorrent, self.sample_dict['raw'], trimmed,
                                            self.cpu, self.parallel, self.size)
                Methods.flag_done(done_trimming)
            else:
                print('Skipping trimming. Already done.')
        else:  # if Illumina reads
            """
            The raw Illumina reads wont typically be demultiplexed by sample and can be cleaned and demultiplexed
            with "process_ragtags" included with Stacks, so trimming with bbduk is disabled.
            """
            print('Processing Illumina reads...')
            if self.read_type == 'se':  # if single-end
                # Trim and demultiplex with stacks
                if not os.path.exists(done_trimming):
                    # It's better no to trim the reads because the de novo process needs reads all the same length
                    # Methods.parallel_trim_reads(Methods.trim_illumina_se, self.sample_dict,
                    #                             cleaned, self.cpu, self.parallel, self.size)
                    Methods.process_radtags_se(self.input, self.barcodes, trimmed, self.enz1, renz_2=self.enz2)
                    Methods.flag_done(done_trimming)
                else:
                    print('Skipping trimming. Already done.')
            else:  # elif self.read_type == 'pe':  # if paired-end
                # Trim and demultiplex
                if not os.path.exists(done_trimming):
                    # It's better no to trim the reads because the de novo process needs reads all the same length
                    # Methods.parallel_trim_reads(Methods.trim_illumina_pe, self.sample_dict,
                    #                             cleaned, self.cpu, self.parallel, self.size)
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

        depth_filt_out = populations + 'populations.depth{}.vcf'.format(self.min_depth)
        maf_filt_out = populations + 'populations.depth{}.maf{}.vcf'.format(
            self.min_depth, self.min_maf)
        missing_filt_out = populations + 'populations.depth{}.maf{}_missing{}.vcf'.format(
            self.min_depth, self.min_maf, self.max_missing)
        ld_stat_out = populations + 'populations.ld_stats.vcf'
        stat_out = populations + 'populations.pop_stats.vcf'
        homo_filt_out = populations + 'populations.depth{}.maf{}_missing{}.homo.vcf'.format(
            self.min_depth, self.min_maf, self.max_missing)

        Methods.vcftools_filter_depth(populations + 'populations.snps.vcf', depth_filt_out, self.min_depth)
        Methods.vcftools_filter_maf(depth_filt_out, maf_filt_out, self.min_maf)
        Methods.vcftools_filter_missing(maf_filt_out, missing_filt_out, self.max_missing)

        # Make tree
        tree = self.out_folder + '/5_tree/'
        Methods.make_folder(tree)

        fastq_out = tree + 'populations.depth{}.maf{}_missing{}.fasta'.format(
            self.min_depth, self.min_maf, self.max_missing)
        fasta_homo_out = tree + 'populations.depth{}.maf{}_missing{}.homo.fasta'.format(
            self.min_depth, self.min_maf, self.max_missing)

        # Only keep homozygous loci
        Methods.filter_out_heterozygous(missing_filt_out, homo_filt_out)

        # Convert VCF to fasta and make tree
        print('Converting VCF to fasta...')
        Methods.vcf2fasta(homo_filt_out, fasta_homo_out)
        Methods.vcf2fasta(missing_filt_out, fastq_out)
        print('Making tree...')
        Methods.make_tree_raxml(fasta_homo_out, tree, self.cpu)
        Methods.make_tree_raxml(fastq_out, tree, self.cpu)

        # Plot trees
        print('Plotting trees...')
        tree_file_in = tree + 'RAxML_bestTree.populations.depth{}.maf{}_missing{}.tree'.format(
            self.min_depth, self.min_maf, self.max_missing)
        tree_render_out = '.'.join(tree_file_in.split('.')[:-1]) + '.png'
        tree_file_homo_in = tree + 'RAxML_bestTree.populations.depth{}.maf{}_missing{}.homo.tree'.format(
            self.min_depth, self.min_maf, self.max_missing)
        tree_render_homo_out = '.'.join(tree_file_homo_in.split('.')[:-1]) + '.png'
        Methods.plot_newick_tree(tree_file_in, tree_render_out)
        Methods.plot_newick_tree(tree_file_homo_in, tree_render_homo_out)

        # Making stats
        print('Making stats...')
        Methods.ld_stat(missing_filt_out, ld_stat_out)
        Methods.vcftools_stats(missing_filt_out, stat_out)

        vcf_dict = Methods.parse_vcf(populations + 'populations.depth{}.maf{}_missing{}.vcf'.format(
            self.min_depth, self.min_maf, self.max_missing))

        # Coverage graph
        print('Making summary stat graphs...')
        stats = self.out_folder + '/6_stats/'
        Methods.make_folder(stats)

        fig_list = Methods.coverage_graph(vcf_dict, stats, self.min_depth, self.min_maf, self.max_missing)
        with open(stats + 'sumstats.html', 'w') as f:
            for fig in fig_list:
                f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


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
                        help='Reference genome for read mapping. '
                             'Mandatory.')
    parser.add_argument('-i', '--input', metavar='/input_folder/',
                        required=True,
                        help='Folder that contains the fastq files. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        help='Folder to hold the result files. '
                             'Mandatory.')
    parser.add_argument('-s', '--size', metavar='64|auto',
                        required=False, default=64,
                        help='Minimum read size to keep after trimming. '
                             'An experimetial auto size selection (use "auto" as argument) is also available. '
                             'It is based on highest peak detection after plotting read length distribution. '
                             'Experimental. Default is 64. Optional.')
    parser.add_argument('-e1', '--enzyme1', metavar='mspI',
                        required=True,
                        help='Restriction enzyme used for the GBS procedure. Mandatory.')
    parser.add_argument('-e2', '--enzyme2', metavar='pstI',
                        required=False,
                        help='Restriction enzyme used for the GBS procedure. Optional.')
    parser.add_argument('-p', '--population-map', metavar='/population_map.tsv',
                        required=True,
                        help='A two-column tab-separated file containing a population map of giving samples.'
                             ' See https://catchenlab.life.illinois.edu/stacks/manual/#popmap for details. '
                             'Mandatory.')
    parser.add_argument('-b', '--barcodes', metavar='/barcodes.tsv',
                        required=True,
                        help='A two-column tab-separated file containing barcode info of giving samples.'
                             ' See https://catchenlab.life.illinois.edu/stacks/manual/#specbc for details. '
                             'Mandatory.')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({}). Optional'.format(max_cpu))
    parser.add_argument('--parallel', metavar='2',
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
    # parser.add_argument('--auto-remove-short-reads',
    #                     required=False,
    #                     action='store_true',  # implies default=False
    #                     help='Auto remove shorter reads based on highest peak detection after plotting read '
    #                          'length distribution. Only works with IonTorrent data. Experimental.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    GBS(arguments)
