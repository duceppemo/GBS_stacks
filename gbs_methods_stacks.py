
# conda install -c bioconda stacks raxml bowtie2 fastp bbmap pysam

import subprocess
import os
import sys
from concurrent import futures
import pathlib
from multiprocessing import cpu_count
from psutil import virtual_memory
import pandas as pd
from collections import defaultdict
from scipy.signal import find_peaks
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np
import gzip


class Methods(object):

    accepted_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
    ion_adapter = 'ATCACCGACTGCCCATAGAGAGG'
    accepted_enz = ['aciI', 'ageI', 'aluI', 'apaLI', 'apeKI', 'apoI', 'aseI', 'bamHI', 'bbvCI', 'bfaI',
                    'bfuCI', 'bgIII', 'bsaHI', 'bspDI', 'bstYI', 'btgI', 'cac8I', 'claI', 'csp6I', 'ddeI',
                    'dpnII', 'eaeI', 'ecoRI', 'ecoRV', 'ecoT22I', 'haeIII', 'hinP1I', 'hindIII', 'hpaII',
                    'hpyCH4IV', 'kpnI', 'mluCI', 'mseI', 'mslI', 'mspI', 'ncoI', 'ndeI', 'ngoMIV', 'nheI',
                    'nlaIII', 'notI', 'nsiI', 'nspI', 'pacI', 'pspXI', 'pstI', 'rsaI', 'sacI', 'sau3AI',
                    'sbfI', 'sexAI', 'sgrAI', 'speI', 'sphI', 'taqI', 'xbaI', 'xhoI']

    @staticmethod
    def check_folder_empty(folder):
        status = 0
        # List content of folder
        dir_content = os.listdir(folder)
        # if folder is not empty and all files have the accepted extensions
        test_file_ext = [x.endswith(tuple(Methods.accepted_extensions)) for x in dir_content]
        if dir_content and any(test_file_ext):
            status = 1
        return status

    @staticmethod
    def check_cpus(requested_cpu, n_proc):
        total_cpu = cpu_count()

        if 1 > requested_cpu > total_cpu:
            requested_cpu = total_cpu
            sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        if 1 > n_proc > total_cpu:
            n_proc = 2
            sys.stderr.write("Number of threads was set to {}".format(2))

        return requested_cpu, n_proc

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def check_mode(ref, mode_ref, mode_denovo):
        if not mode_ref and not mode_denovo:
            raise Exception('Please choose your analysis mode: "--referenced" or "--de-novo".')
        if mode_ref:
            if not ref:
                raise Exception('Please specify a reference fasta file for your referenced analysis.')

    @staticmethod
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder):
        sample_dict = defaultdict(list)
        # sample_dict = dict()
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(Methods.accepted_extensions)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    sample = filename.split('.')[0].split('_R1')[0].split('_R2')[0]
                    if '_R1' in filename:
                        sample_dict[sample].insert(0, file_path)
                    elif '_R2' in filename:
                        sample_dict[sample].insert(1, file_path)
                    else:
                        sample_dict[sample].insert(0, file_path)
        if not sample_dict:
            raise Exception('Sample dictionary empty!')

        return sample_dict

    @staticmethod
    def check_map(sample_dict, map_file):
        map_file_sample_list = list()
        with open(map_file, 'r') as f:
            for line in f:
                map_file_sample_list.append(line.split('\t')[0])  # First column is sample name

        # Check if sample from map file have sequencing data
        missing_list = list()
        for sample, info in sample_dict.items():
            if sample not in map_file_sample_list:
                missing_list.append(sample)
        if missing_list:
            raise Exception('The following samples from your population map file are missing corresponding'
                            ' sequencing files:{}'.format(', '.join(missing_list)))
        # Check if some fastq are not in the population map file
        # Actually, I think this is OK. I will just waste time processing them, but should work with Stacks.

    @staticmethod
    def check_enzymes(enz_list):
        for enz in enz_list:
            if enz not in Methods.accepted_enz:
                raise Exception('Only the following restriction enzymes are accepted: {}'
                                .format(', '.join(Methods.accepted_enz)))

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def flag_done(flag_file):
        with open(flag_file, 'w') as f:
            pass

    @staticmethod
    def trim_illumina_se(r1, output_folder, cpu):
        sample = os.path.basename(r1).split('_R1')[0]

        Methods.make_folder(output_folder)

        cmd = ['fastp',
               '--thread', str(cpu),
               '--in1', r1,
               '--out1', output_folder + '/' + sample + 'fastq.gz',
               '--length_required', str(64),
               '--cut_right',
               '--html', output_folder + '/' + sample + '.html']

        print('\t{}'.format(sample))
        subprocess.run(cmd, stderr=subprocess.DEVNULL)

    @staticmethod
    def trim_illumina_pe(r1, output_folder, cpu):
        r2 = r1.replace('_R1', '_R2')
        sample = os.path.basename(r1).split('_R1')[0]

        Methods.make_folder(output_folder)

        cmd = ['fastp',
               '--thread', str(cpu),
               '--in1', r1,
               '--in2', r2,
               '--out1', output_folder + '/' + sample + '_R1.fastq.gz',
               '--out2', output_folder + '/' + sample + '_R2.fastq.gz',
               '--length_required', str(64),
               '--cut_right',
               '--html', output_folder + '/' + sample + '.html']

        print('\t{}'.format(sample))
        subprocess.run(cmd, stderr=subprocess.DEVNULL)

    @staticmethod
    def trim_iontorrent(r1, trimmed_folder, cpu):
        sample = os.path.basename(r1).split('_R1')[0]
        sample = sample.split('.')[0]

        cmd = ['bbduk.sh',
               'in={}'.format(r1),
               'out={}'.format(trimmed_folder + sample + '.fastq.gz'),
               'literal={}'.format(Methods.ion_adapter),
               'k=21', 'mink=15', 'hdist=1', 'ktrim=r', 'trimq=6', 'minlength=64',
               'overwrite=t', 'threads={}'.format(cpu)]

        print('\t{}'.format(sample))
        # subprocess.run(cmd, stderr=subprocess.DEVNULL)
        log_file = trimmed_folder + sample + '_bbduk.log'
        with open(log_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

    @staticmethod
    def read_length_dist(r1):
        """
        Make read length start of all reads combined, plot the distribution and find the peak.
        Returns the optimal read length for trimming the reads for the ustacks.
        """
        sample = os.path.basename(r1).split('_R1')[0]
        sample = sample.split('.')[0]

        len_dict = dict()

        cmd = ['readlength.sh',
               'in={}'.format(r1),
               'bin=5']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        len_dist = p.communicate()[0].decode('utf-8')
        lines = len_dist.split('\n')
        for line in lines:
            if line.startswith('#'):
                continue
            elif not line:
                continue
            else:
                fields = line.split('\t')[:2]
                len_dict[fields[0]] = fields[1]

        return len_dict

    @staticmethod
    def parallel_read_length_dist(sample_dict, output_folder, cpu):
        master_len_dict = dict()
        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = (path_list[0] for sample, path_list in sample_dict.items())
            for results in executor.map(Methods.read_length_dist,
                                        [path_list[0] for sample, path_list in sample_dict.items()]):
                for k, v in results.items():
                    if k in master_len_dict:
                        master_len_dict[int(k)] += int(v)
                    else:
                        master_len_dict[int(k)] = int(v)
        df = pd.DataFrame.from_dict(master_len_dict, orient='index')  # , columns=['Length', 'Count'])
        df.index.names = ['Length']
        df.rename(columns={df.columns[0]: 'Count'}, inplace=True)

        return df

    @staticmethod
    def find_peak(df, output_folder):
        # plot the distribution
        x = df.loc[:, 'Count'].values
        y = df.index.values
        xs = np.sort(x)
        ys = np.array(y)[np.argsort(x)]

        peaks, properties = find_peaks(x, prominence=1)
        highest_peak = int(properties['prominences'].max())
        p = np.interp(highest_peak, ys, xs)
        size_to_keep = int(df[df['Count'] == int(p)].index.values)
        fig, ax = plt.subplots()
        g = plt.plot(x)
        plt.plot(peaks, x[peaks], "x")
        plt.tight_layout()
        fig.savefig(output_folder + "/peaks.png")
        plt.close()

        return size_to_keep

    @staticmethod
    def trim_iontorrent_size(r1, trimmed_folder, cpu, size):
        """
        trim all reads to 100bp
        """
        sample = os.path.basename(r1).split('_R1')[0]
        sample = sample.split('.')[0]

        cmd = ['bbduk.sh',
               'in={}'.format(r1),
               'out={}'.format(trimmed_folder + sample + '.fastq.gz'),
               'literal={}'.format(Methods.ion_adapter),
               'k=21', 'mink=15', 'hdist=1', 'ktrim=r', 'trimq=6',
               'minlength={}'.format(size), 'forcetrimright={}'.format(size-1),
               'overwrite=t', 'threads={}'.format(cpu)]

        print('\t{}'.format(sample))
        # subprocess.run(cmd, stderr=subprocess.DEVNULL)
        log_file = trimmed_folder + sample + '_bbduk.log'
        with open(log_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_trim_reads(trim_function, sample_dict, output_folder, cpu, parallel, **kwargs):
        print('Trimming reads...')
        Methods.make_folder(output_folder)
        test = kwargs['size']

        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            if kwargs:
                args = ((path_list[0], output_folder, int(cpu / parallel), kwargs['size'])
                        for sample, path_list in sample_dict.items())
            else:
                args = ((path_list[0], output_folder, int(cpu / parallel))
                        for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: trim_function(*x), args):
                pass

    @staticmethod
    def index_bowtie2(ref, prefix, cpu):
        print('Indexing reference genome...')
        cmd = ['bowtie2-build', '--threads', str(cpu), ref, prefix]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_bowtie2_se(ref, fastq_file, cpu, output_bam):
        ref_index = '.'.join(ref.split('.')[:-1])
        sample = '.'.join(os.path.basename(fastq_file).split('.')[:-1])
        if fastq_file.endswith('.gz'):
            sample = '.'.join(sample.split('.')[:-1])

        print('\t{}'.format(sample))
        bowtie2_align_cmd = ['bowtie2',
                             '-x', ref_index,
                             '-U', fastq_file,
                             '--threads', str(cpu),
                             '--rg-id', sample,
                             '--rg', 'SM:{}'.format(sample)]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4', '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index', '-c',
                              output_bam]

        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.communicate()
        # TODO: log bowtie2 output (stats) to file

        # index cram file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_map_bowtie2_se(output_folder, ref, sample_dict, cpu, parallel):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((ref, path_list[0], int(cpu/parallel), output_folder + '/' + sample + '.bam')
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2_se(*x), args):
                pass

    @staticmethod
    def map_bowtie2_pe(ref, r1, r2, cpu, output_bam):
        ref_index = '.'.join(ref.split('.')[:-1])
        sample = '.'.join(os.path.basename(r1).split('.')[:-1])
        if r1.endswith('.gz'):
            sample = '.'.join(sample.split('.')[:-1])

        print('\t{}'.format(sample))
        bowtie2_align_cmd = ['bowtie2',
                             '--xeq ',
                             '-x', ref_index,
                             '-1', r1,
                             '-2', r2,
                             '--threads', str(cpu),
                             '--rg-id', sample,
                             '--rg', 'SM:{}'.format(sample)]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-F', '4',
                             '-h',
                             '-T', ref,
                             '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        # samtools can only index chromosomes up to 512M bp.
        samtools_index_cmd = ['samtools', 'index', '-c',
                              output_bam]

        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2.stdout.close()
        p3.communicate()
        # TODO: log bowtie2 output (stats) to file

        # index cram file
        subprocess.run(samtools_index_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def parallel_map_bowtie2_pe(output_folder, ref, sample_dict, cpu, parallel):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((ref, path_list[0], path_list[1], int(cpu/parallel), output_folder + '/' + sample + '.bam')
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2_pe(*x), args):
                pass

    @staticmethod
    def process_radtags_se(input_fastq_folder, barcode_tsv, output_folder, enz, **kwargs):
        Methods.make_folder(output_folder)
        extra_arg_list = ['--renz_2']
        if kwargs:
            for k, v in kwargs.items():
                extra_arg_list.append(v)
        cmd = ['process_radtags',
               '-p', input_fastq_folder,
               '-b', barcode_tsv,
               '-o', output_folder,
               '--renz_1', enz,
               '--clean', '--quality', '--rescue'] + extra_arg_list
        subprocess.run(cmd)

    @staticmethod
    def process_radtags_pe(input_fastq_folder, barcode_tsv, output_folder, enz, **kwargs):
        Methods.make_folder(output_folder)
        cmd = ['process_radtags',
               '-p', input_fastq_folder,
               '--paired',
               '-b', barcode_tsv,
               '-o', output_folder,
               '--renz_1', enz,
               kwargs,
               '--clean', '--quality', '--rescue']
        subprocess.run(cmd)

    @staticmethod
    def ustacks(r1, output_folder, unique_id, cpu):
        Methods.make_folder(output_folder)
        cmd = ['ustacks',
               '-f', r1,
               '-i', unique_id,
               '-o', output_folder,
               '-p', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def parallel_ustacks(sample_dict, output_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((path_list[0], output_folder, str(i), int(cpu / parallel))
                    for i, (sample, path_list) in enumerate(sample_dict.items()))
            for results in executor.map(lambda x: Methods.ustacks(*x), args):
                pass

    @staticmethod
    def cstacks(input_file, stacks_folder, cpu):
        sample = stacks_folder + os.path.basename(input_file).split('.')[0]
        cmd = ['cstacks',
               '-s', sample,
               '-o', stacks_folder,
               '-p', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def parallel_cstacks(sample_dict, output_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((path_list[0], output_folder, int(cpu/parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.cstacks(*x), args):
                pass

    @staticmethod
    def sstacks(input_file, stacks_folder, cpu):
        sample = stacks_folder + os.path.basename(input_file).split('.')[0]
        cmd = ['sstacks',
               '-c', stacks_folder,
               '-s', sample,
               '-o', stacks_folder,
               '-p', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def parallel_sstacks(sample_dict, output_folder, cpu, parallel):
        with futures.ThreadPoolExecutor(max_workers=parallel) as executor:
            args = ((path_list[0], output_folder, int(cpu/parallel))
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.sstacks(*x), args):
                pass

    @staticmethod
    def tsv2bam_se(stacks_folder, pop_map, cpu):
        cmd = ['tsv2bam',
               '-P', stacks_folder,
               '-M', pop_map,
               '-t', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def tsv2bam_pe(stacks_folder, pop_map, fastq_folder, cpu):
        cmd = ['tsv2bam',
               '-P', stacks_folder,
               '-M', pop_map,
               '-R', fastq_folder,
               '-t', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def call_snps_gstacks(bam_folder, pop_map, output_folder, cpu):
        cmd = ['gstacks',
               '-I', bam_folder,
               '-M', pop_map,
               '-O', output_folder,
               '-t', str(cpu)]
        subprocess.run(cmd)

    @staticmethod
    def make_pop_stats(gstacks_folder, output_folder, pop_map, cpu):
        Methods.make_folder(output_folder)

        # num_lines = sum(1 for line in open(pop_map, 'r') if line.rstrip() != '')
        cmd = ['populations',
               '--in-path', gstacks_folder,
               '--out-path', output_folder,
               '--popmap', pop_map,
               '--threads', str(cpu),
               '--vcf',
               '--write-single-snp']
        subprocess.run(cmd)

    @staticmethod
    def vcftools_filter_snp(input_vcf, output_vcf):
        # Only keep biallelic SNP loci
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--remove-indels',
               '--remove-filtered-all',
               '--min-alleles', str(2), '--max-alleles', str(2),
               '--recode']
        subprocess.run(cmd)

        # Remove "recode" from the output file name
        os.replace(output_vcf + '.recode.vcf', output_vcf)

    # @staticmethod
    # def vcftools_filter_qual(input_vcf, output_vcf):
    #     cmd = ['vcftools',
    #            '--vcf', input_vcf,
    #            '--out', output_vcf,
    #            '--minQ', str(20),
    #            '--recode']
    #     subprocess.run(cmd)
    #
    #     # Remove "recode" from the output file name
    #     os.replace(output_vcf + '.recode.vcf', output_vcf)

    @staticmethod
    def vcftools_filter_maf(input_vcf, output_vcf, min_maf):
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--maf', str(min_maf),
               '--recode']
        subprocess.run(cmd)

        # Remove "recode" from the output file name
        os.replace(output_vcf + '.recode.vcf', output_vcf)

    @staticmethod
    def vcftools_filter_depth(input_vcf, output_vcf, min_depth):
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--min-meanDP', str(min_depth),
               '--recode']
        subprocess.run(cmd)

        # Remove "recode" from the output file name
        os.replace(output_vcf + '.recode.vcf', output_vcf)

    @staticmethod
    def vcftools_filter_missing(input_vcf, output_vcf, max_missing):
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--max-missing', str(max_missing),
               '--recode']
        subprocess.run(cmd)

        # Remove "recode" from the output file name
        os.replace(output_vcf + '.recode.vcf', output_vcf)

    @staticmethod
    def vcftools_stats(input_vcf, output_vcf):
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--site-mean-depth']
        subprocess.run(cmd)

    @staticmethod
    def ld_filering(input_vcf, output_vcf):
        # Output sites with Linkage desequilibrium
        cmd = ['vcftools',
               '--vcf', input_vcf,
               '--out', output_vcf,
               '--hap-r2', '--min-r2', str(0.5)]

        ld_file = output_vcf + '.hap.ld'
        subprocess.run(cmd)

    @staticmethod
    def filter_out_heterozygous(input_vcf, output_vcf):
        with open(output_vcf, 'w') as out_fh:
            with open(input_vcf, 'r') as in_fh:
                match = ['0/1', '1/0']
                for line in in_fh:
                    if any(x in line for x in match):
                        continue
                    else:
                        out_fh.write(line)

    @staticmethod
    def vcf2fasta(input_vcf, output_fasta):
        iupac_dict = {'AG': 'R',
                      'CT': 'Y',
                      'GC': 'S',
                      'AT': 'W',
                      'GT': 'K',
                      'AC': 'M'}

        with open(output_fasta, 'w') as out_fh:
            with open(input_vcf, 'r') as in_fh:
                for line in in_fh:
                    if line.startswith('##'):
                        # out_fh.write(line)
                        continue
                    elif line.startswith('#CHROM'):
                        sample_list = line.rstrip().split('\t', 9)[9].split('\t')
                        df = pd.DataFrame(columns=sample_list)
                    else:
                        # Split data lines into 10 fields, the last one is the samples info
                        field_list = line.rstrip().split('\t', 9)
                        # Dictionary to convert 0/0 and 1/1 geno to REF or ALT call
                        ref = field_list[3]
                        alt = field_list[4]
                        try:
                            iupac = iupac_dict[ref + alt]
                        except KeyError:
                            iupac = iupac_dict[alt + ref]

                        geno_dict = {'0/0': ref,
                                     '1/1': alt,
                                     '0/1': iupac,
                                     '1/0': iupac,
                                     './.': '-'}
                        # Split the last field (sample info) and only keep the genotype
                        df.loc[len(df)] = [geno_dict[x.split(':')[0]] for x in field_list[9].split('\t')]

                # Convert dataframe to fasta output
                fasta_dict = df.to_dict(orient='list')
                for sample, seq_list in fasta_dict.items():
                    out_fh.write('>{}\n{}\n'.format(sample, ''.join(seq_list)))

    @staticmethod
    def make_tree_raxml(aligned_fasta, output_folder, cpu):
        Methods.make_folder(output_folder)

        cmd = ['raxmlHPC-PTHREADS-AVX2',
               '-s', aligned_fasta,
               '-w', output_folder,
               '-n', '{}.tree'.format('.'.join(os.path.basename(aligned_fasta).split('.')[:-1])),
               '-m', 'GTRCAT',
               '-N', str(100),
               '-d', '-f', 'a', '-T', str(cpu),
               '-x', str(1234), '-p', str(123)]

        subprocess.run(cmd)

    @staticmethod
    def parse_vcf(vcf_file):
        vcf_dict = defaultdict()
        with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    sample_list = line.rstrip().split('\t', 9)[9].split('\t')
                    # Convert list to dictionary
                    # vcf_dict = {i: '' for i in sample_list}
                else:
                    # Split data lines into 10 fields, the last one is the samples info
                    field_list = line.rstrip().split('\t', 9)
                    # Dictionary to convert 0/0 and 1/1 geno to REF or ALT call
                    ref = field_list[3]
                    alt = field_list[4]

                    # Split the last field (sample info) and only keep the genotype
                    for i, sample_info in enumerate(field_list[9].split('\t')):
                        geno_fieds = sample_info.split(':')
                        GT = geno_fieds[0]
                        DP = geno_fieds[1]

                        if sample_list[i] not in vcf_dict:
                            vcf_dict[sample_list[i]] = {'GT': [], 'DP': []}
                        vcf_dict[sample_list[i]]['GT'].append(GT)
                        vcf_dict[sample_list[i]]['DP'].append(DP)
        return vcf_dict

    @staticmethod
    def coverage_graph(vcf_dict, output_folder):
        # graph_dict = dict()
        # n_loci = 0
        # for sample, info_dict in vcf_dict.items():
        #     if sample not in graph_dict:
        #         graph_dict[sample] = {'cov': {}, 'miss': {}, 'hetero': {}}
        #
        #     c = 0
        #     for x in info_dict['DP']:
        #         n_loci += 1
        #         try:
        #             c += int(x)
        #         except ValueError:
        #             continue
        #     graph_dict[sample]['cov'] = c
        #
        #     hetero = 0
        #     missing = 0
        #     for x in info_dict['GT']:
        #         if x == './.':
        #             missing += 1
        #         elif x == ('0/1' or '1/0'):
        #             hetero += 1
        #     graph_dict[sample]['miss'] = missing
        #     graph_dict[sample]['hetero'] = hetero

        # Coverage per sample
        df = pd.DataFrame(columns=['Sample', 'GT', 'DP'], index=None)
        for sample, info_dict in vcf_dict.items():
            for i in range(len(info_dict['GT'])):
                gt = info_dict['GT'][i]
                dp = info_dict['DP'][i]

                df = df.append({'Sample': sample, 'GT': gt, 'DP': dp}, ignore_index=True)

        df.replace(r'\.|\./\.', np.nan, regex=True, inplace=True)
        df1 = df.loc[:, ['Sample', 'DP']].dropna()
        df1['DP'] = df1['DP'].astype(int)
        fig = px.violin(df1, x='Sample', y='DP')
        fig.write_html(output_folder + '/' + 'coverage.html', auto_open=False)

        # % missing per sample
        n_loci = len(df)
        df2 = df.GT.isnull().groupby(df['Sample']).sum().astype(int).reset_index(name='% Missing')
        df2['% Missing'] = df2['% Missing'].apply(lambda x: int(x / n_loci * 100))
        fig = px.bar(df2, x='Sample', y='% Missing')
        fig.write_html(output_folder + '/' + 'missing.html', auto_open=False)

        # % Hetorozygote per sample
        df['GT'] = df['GT'].astype('category')
        df3 = df.groupby(['Sample', 'GT']).size().reset_index(name="% Hetero")
        df3['% Hetero'] = df3['% Hetero'].apply(lambda x: int(x / n_loci * 100))
        fig = px.bar(df3, x='Sample', y='% Hetero')
        fig.write_html(output_folder + '/' + 'hetero.html', auto_open=False)

        # Markers LD distribution???

        pass
