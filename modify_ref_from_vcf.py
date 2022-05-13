#!/usr/local/env python3

import gzip
from argparse import ArgumentParser
from itertools import groupby


__author__ = 'duceppemo'
__version__ = 'v0.1'


class RefExtractor(object):
    def __init__(self, args):
        # Arguments
        self.ref = args.reference
        self.vcf = args.variant
        self.output = args.output_fasta

        # Run
        vcf_dict = Methods.parse_vcf(self.vcf)
        Methods.extract_regions_iter(self.ref, vcf_dict, self.output)


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
                ref = field_list[3]
                alt = field_list[4]
                if chrom not in vcf_dict:
                    vcf_dict[chrom] = [(pos, ref, alt)]
                else:
                    if pos not in vcf_dict[chrom]:  # Avoid duplicated positions
                        vcf_dict[chrom].append((pos, ref, alt))

        return vcf_dict

    @staticmethod
    def extract_regions_iter(ref, vcf_dict, out):
        counter = 1

        # https://www.biostars.org/p/710/
        with open(out, 'w') as fh_out:
            print('Parsing reference fasta file...')
            with gzip.open(ref, 'rt') if ref.endswith('.gz') else open(ref, 'r') as fh_in:
                # Create iterator
                faiter = (x[1] for x in groupby(fh_in, lambda line: line[0] == '>'))

                for header in faiter:
                    # Writing form header
                    fasta_header = header.__next__()
                    fh_out.write(fasta_header)
                    chrom = fasta_header[1:].rstrip().split()[0]  # drop the ">"

                    # join all sequence lines to one.
                    full_seq = ''.join(s.rstrip() for s in faiter.__next__())

                    if chrom in vcf_dict:  # contig (or "chromosome" needs to be changed)
                        for pos_tuple in vcf_dict[chrom]:
                            pos, ref, alt = pos_tuple
                            index = int(pos) - 1  # position 1 is index 0

                            # Change the reference
                            ref_og = full_seq[index]
                            full_seq = full_seq[:index] + ref + full_seq[index + 1:]
                            ref_new = full_seq[index]
                            t = 1

                    # Write sequence to file.
                    fh_out.write('{}\n'.format(full_seq))


if __name__ == "__main__":
    parser = ArgumentParser(description='Modify a fasta file to match the REF call from a VCF file.')
    parser.add_argument('--reference', '-ref', metavar='reference.fasta',
                        required=True,
                        type=str,
                        help='Reference fasta file. Must be the same that was used to create the VCF file. Mandatory.')
    parser.add_argument('--variant', '-vcf', metavar='variants.vcf',
                        required=True,
                        type=str,
                        help='VCF file for which reference regions will be extracted for each position. Mandatory.')
    parser.add_argument('--output-fasta', '-out', metavar='modified_ref.fasta',
                        type=str,
                        required=True,
                        help='Fasta output file containing the regions. Mandatory')
    # Get the arguments into an object
    arguments = parser.parse_args()
    RefExtractor(arguments)
