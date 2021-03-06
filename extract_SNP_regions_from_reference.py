#!/usr/local/env python3

import gzip
from argparse import ArgumentParser
from itertools import groupby


__author__ = 'duceppemo'
__version__ = 'v0.1'


class Ref_Extractor(object):
    def __init__(self, args):
        # Arguments
        self.ref = args.reference
        self.vcf = args.variant
        self.bp = args.base_pairs
        self.output = args.output_fasta

        # Run
        vcf_dict = Methods.parse_vcf(self.vcf)
        Methods.extract_regions_iter(self.ref, vcf_dict, self.bp, self.output)


class Methods(object):
    @staticmethod
    def parse_vcf(vcf_file):
        print('Parsing VCF file...')
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
                    if pos not in vcf_dict[chrom]:  # Avoid duplicated positions
                        vcf_dict[chrom].append(pos)

        return vcf_dict

    @staticmethod
    def extract_regions_iter(ref, vcf_dict, bp, out):
        # https://www.biostars.org/p/710/
        with open(out, 'w') as fh_out:
            print('Parsing reference fasta file...')
            with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as fh_in:
                # Create iterator
                faiter = (x[1] for x in groupby(fh_in, lambda line: line[0] == '>'))

                for header in faiter:
                    # drop the ">"
                    chrom = header.__next__()[1:].rstrip().split()[0]
                    # join all sequence lines to one.
                    full_seq = ''.join(s.rstrip() for s in faiter.__next__())
                    if chrom in vcf_dict:
                        for pos in vcf_dict[chrom]:
                            index = int(pos) - 1  # position 1 is index 0
                            start = index - int(bp)
                            if start <= 0:  # if the SNP is closer to beginning of sequence than the range required
                                start = 0
                            stop = index + int(bp)
                            if stop > len(full_seq):  # if the SNP is closer to end of sequence than the range required
                                stop = len(full_seq)

                            header = '>{}:{}-{} SNP at {}'.format(chrom, start + 1, stop + 1, index + 1)  # index -> position
                            extracted_region = full_seq[start:stop + 1]
                            fh_out.write('{}\n{}\n'.format(header, extracted_region))

                    # yield tuple(header, seq)


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
                        required=False,
                        type=int,
                        default=200,
                        help='Number of base pairs to extract on each side of the SNP. Defaults is 200. Mandatory.')
    parser.add_argument('--output-fasta', '-out', metavar='SNP_regions.fasta',
                        type=str,
                        required=True,
                        help='Fasta output file containing the regions. Mandatory')
    # Get the arguments into an object
    arguments = parser.parse_args()
    Ref_Extractor(arguments)
