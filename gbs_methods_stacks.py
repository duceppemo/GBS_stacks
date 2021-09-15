
# conda install -c bioconda stacks raxml bowtie2 fastp bbmap pysam

import subprocess
import os
import sys
import pysam
from concurrent import futures
import pathlib
from multiprocessing import cpu_count
from psutil import virtual_memory
import pandas as pd


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
    def check_cpus(requested_cpu):
        total_cpu = cpu_count()

        if requested_cpu:
            if requested_cpu > total_cpu:
                requested_cpu = total_cpu
                sys.stderr.write("Number of threads was set higher than available CPUs ({})".format(total_cpu))
                sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        else:
            requested_cpu = total_cpu

        return requested_cpu

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
    def make_folder(folder):
        # Will create parent directories if don't exist and will not return error if already exists
        pathlib.Path(folder).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def get_files(in_folder, sample_dict):
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
    def trim_illumina_se(r1, trimmed_folder, cpu):
        sample = os.path.basename(r1).split('_R1')[0]

        cmd = ['fastp',
               '--thread', str(cpu),
               '--in1', r1,
               '--out1', trimmed_folder + '/' + sample + 'fastq.gz',
               '--length_required', str(64),
               '--cut_right',
               '--html', trimmed_folder + '/' + sample + '.html']

        print('\t{}'.format(sample))
        subprocess.run(cmd, stderr=subprocess.DEVNULL)

    @staticmethod
    def trim_illumina_pe(r1, trimmed_folder, cpu):
        r2 = r1.replace('_R1', '_R2')
        sample = os.path.basename(r1).split('_R1')[0]

        cmd = ['fastp',
               '--thread', str(cpu),
               '--in1', r1,
               '--in2', r2,
               '--out1', trimmed_folder + '/' + sample + '_R1.fastq.gz',
               '--out2', trimmed_folder + '/' + sample + '_R2.fastq.gz',
               '--length_required', str(64),
               '--cut_right',
               '--html', trimmed_folder + '/' + sample + '.html']

        print('\t{}'.format(sample))
        subprocess.run(cmd, stderr=subprocess.DEVNULL)

    @staticmethod
    def trim_iontorrent(r1, trimmed_folder, cpu):
        sample = os.path.basename(r1).split('_R1')[0]

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
    def parallel_trim_reads(trim_function, sample_dict, output_folder, cpu):
        print('Trimming reads...')
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            args = ((path_list[0], output_folder, int(cpu/4))
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
    def parallel_map_bowtie2_se(output_folder, ref, sample_dict, cpu):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            args = ((ref, path_list[0], int(cpu/4), output_folder + '/' + sample + '.bam')
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
    def parallel_map_bowtie2_pe(output_folder, ref, sample_dict, cpu):
        # Index reference genome if not already done
        ref_index = '.'.join(ref.split('.')[:-1])
        if not os.path.exists(ref_index + '.1.bt2') and not os.path.exists(ref_index + '.1.bt2l'):
            Methods.index_bowtie2(ref, ref_index, cpu)

        print('Mapping reads...')
        Methods.make_folder(output_folder)

        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            args = ((ref, path_list[0], path_list[1], int(cpu/4), output_folder + '/' + sample + '.bam')
                    for sample, path_list in sample_dict.items())
            for results in executor.map(lambda x: Methods.map_bowtie2_pe(*x), args):
                pass

    @staticmethod
    def process_radtags_se(input_fastq_folder, barcode_tsv, output_folder, enz, **kwargs):
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
    def call_snps_gstacks(cram_folder, pop_map, output_folder, cpu):
        cmd = ['gstacks',
               '-I', cram_folder,
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
               '--min-samples-overall', str(1.0),
               '--plink',
               # '--fstats', '--vcf', '--structure', '--fasta-loci', '--plink',
               '--min-maf', str(0.05),
               '--write-single-snp']

        subprocess.run(cmd)

    @staticmethod
    def filter_vcf(input_vcf, output_vcf):
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
        with open(output_fasta, 'w') as out_fh:
            with open(input_vcf, 'r') as in_fh:
                vcf_dict = dict()
                for line in in_fh:
                    if line.startswith('##'):
                        # out_fh.write(line)
                        continue
                    elif line.startswith('#CHROM'):
                        sample_list = line.rstrip().split('\t', 9)[9].split('\t')
                        df = pd.DataFrame(columns=sample_list)
                        # for s in sample_list:
                        #     vcf_dict[s] = []

                    else:
                        # Split data lines into 10 fields, the last one is the samples info
                        field_list = line.rstrip().split('\t', 9)
                        # Dictionary to convert 0/0 and 1/1 geno to REF or ALT call
                        geno_dict = {'0/0': field_list[3],
                                     '1/1': field_list[4],
                                     '0/1': '.',
                                     '1/0': '.'}
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