#!/usr/local/env python3

import os
import gzip
from argparse import ArgumentParser
from collections import defaultdict


__author__ = 'duceppemo'
__version__ = 'v0.1'


class Ref_Extractor(object):
    def __init__(self, args):
        # Arguments
        self.ref = args.reference
        self.vcf = args.variant
        self.bp = args.base_pairs
        self.output = args.output

        # Run
        vcf_dict = Methods.parse_vcf(self.vcf)
        Methods.extract_regions(self.ref, vcf_dict, self.bp, self.output)


class Methods(object):
    @staticmethod
    def parse_vcf(vcf_file):
        vcf_dict = dict()
        with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                if not line:
                    continue
                field_list = line.split('\t')
                chrom = field_list[0]
                pos = field_list[1]
                if chrom not in vcf_dict:
                    vcf_dict[chrom] = [pos]
                else:
                    vcf_dict[chrom].append(pos)

        return vcf_dict

    @staticmethod
    def extract_regions(ref, vcf_dict, bp, out):
        with open(out, 'w') as fh_out:
            ref_dict = dict()
            with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as fh_in:
                for line in fh_in:
                    line = line.rstrip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        chrom = line.split()[0].replace('>', '')
                        ref_dict[chrom] = ''
                    else:
                        ref_dict[chrom] += line
            for chrom, pos in vcf_dict.items():
                index = int(pos) - 1
                start = index - int(bp)
                if start <= 0:
                    start = 0
                stop = index + int(bp)
                if stop > len(ref_dict[chrom]):
                    stop = len(ref_dict[chrom])

                header = '>{}:{}-{}'.format(chrom, start, stop)
                seq = ref_dict[chrom][start:stop]
                fh_out.write('{}\n{}\n'.format(header, seq))


if __name__ == "__main__":
    parser = ArgumentParser(description='Extract flanking regions of SNPs from a VCF file using the corresponding '
                                        'reference fasta file.')
    parser.add_argument('--reference', '-ref', metavar='reference.fasta',
                        required=True,
                        type=str,
                        help='Reference fasta file. Must be the same that was used to create the VCF file. Mandatory.')
    parser.add_argument('--variant', '-vcf', metavar='variants.vcf',
                        required=True,
                        type=str,
                        help='VCF file for which reference regions will be extracted for each position. Mandatory.')
    parser.add_argument('--base-pairs', '-bp', metavar='200',
                        required=True,
                        type=int,
                        default=200,
                        help='Number of base pairs to extract on each side of the SNP. Defaults is 200. Mandatory.')
    parser.add_argument('--output-fasta', '-out', metavar='SNP_regions.fasta',
                        type=str,
                        required=True,
                        help='Fasta output file contraining the regions. Mandatory')
    # Get the arguments into an object
    arguments = parser.parse_args()
    Ref_Extractor(arguments)
