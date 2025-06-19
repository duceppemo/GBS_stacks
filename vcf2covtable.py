from typing import Any
import re
from io import TextIOBase
# import plotly.graph_objects as go
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import pandas as pd
from collections import defaultdict
import numpy as np
import os
import sys
from concurrent import futures

from fontTools.misc.cython import returns

pd.set_option('future.no_silent_downcasting', True)


if len(sys.argv) != 2:
    print('Usage: pythong vcf_depth_table.py <samples.vcf>')
    sys.exit()

my_vcf = sys.argv[1]

# def parse_vcf(vcf_file):
#     vcf_dict = defaultdict()
#     with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
#         for line in f:
#             if line.startswith('##'):
#                 continue
#             elif line.startswith('#CHROM'):
#                 sample_list = line.rstrip().split('\t', 9)[9].split('\t')
#                 # Convert list to dictionary
#                 # vcf_dict = {i: '' for i in sample_list}
#             else:
#                 # Split data lines into 10 fields, the last one is the samples info
#                 field_list = line.rstrip().split('\t', 9)
#                 # Dictionary to convert 0/0 and 1/1 geno to REF or ALT call
#                 ref = field_list[3]
#                 alt = field_list[4]
#
#                 # Split the last field (sample info) and only keep the genotype
#                 for i, sample_info in enumerate(field_list[9].split('\t')):
#                     # GT:DP:AD:GQ:GL
#                     geno_fieds = sample_info.split(':')
#                     if len(geno_fieds) == 1:  # ggeno_fieds[0] == './.'
#                         DP = 0
#                     else:
#                         DP = geno_fieds[1]
#
#                     if sample_list[i] not in vcf_dict:
#                         vcf_dict[sample_list[i]] = {'DP': []}
#                     vcf_dict[sample_list[i]]['DP'].append(DP)
#     return vcf_dict

def parse_vcf(my_vcf_file, my_sample_list):
    my_dict = dict()
    with gzip.open(my_vcf_file, 'rt') if my_vcf_file.endswith('.gz') else open(my_vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                # Split data lines into 10 fields, the last one is the samples info
                field_list = line.rstrip().split('\t', 9)
                chrom = field_list[0]
                pos = field_list[1]
                coord = '{}:{}'.format(chrom, pos)

                # Split the last field (sample info) and only keep the genotype
                for i, sample_info in enumerate(field_list[9].split('\t')):
                    try:
                        # GT:DP:AD:GQ:GL
                        dp = int(sample_info.split(':')[1])
                    except (IndexError, ValueError):
                        dp = 0
                    if my_sample_list[i] not in my_dict:
                        my_dict[my_sample_list[i]] = dict()
                    my_dict[my_sample_list[i]][coord] = dp
        return my_dict

def get_sample_list(vcf_file):
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        my_sample_list: list[Any] = list()
        for line in f:
            if line.startswith('#CHROM'):
                 my_sample_list = line.rstrip().split('\t', 9)[9].split('\t')
            else:
                continue
        return my_sample_list

def process_line(line, my_sample_list):
    """Function to process a single line."""
    my_dict: dict[str, list[int]] = dict()

    if line.startswith('#'):
        pass
    else:
        # Split data lines into 10 fields, the last one is the samples info
        field_list = line.rstrip().split('\t', 9)

        # Split the last field (sample info) and only keep the genotype
        for i, sample_info in enumerate(field_list[9].split('\t')):
            try:
                # GT:DP:AD:GQ:GL
                dp = int(sample_info.split(':')[1])
            except IndexError:
                dp = 0
            if my_sample_list[i] not in my_dict:
                my_dict[my_sample_list[i]] = list()
            my_dict[my_sample_list[i]].append(dp)
    return my_dict

def chunk_read(f_obj: TextIOBase, sentinel: str, max_sentinel: int):
    """Read a file object in chunks
    Read the file object line by line. Each time a sentinel is detected, we increment
    a count. Once the count reaches max_sentinel, we have gatherered the required
    chunk and yield it.
    The function is inspired by this SO answer:
    https://stackoverflow.com/a/42964612/9723036
    NOTE: during chunking, we remove all the white spaces and tabs to reduce the
    memory load.
    :param f_obj: A file object from opening a text file.
    :type f_obj: TextIOBase
    :param sentinel: A string pattern (regex supported) to recognize a specific
        line.
    :type sentinel: str
    :param max_sentinel: Max number of appearance of sentinels allowed in a chunk.
        This is equivalent to a chunk size, but more meaningful than based on only
        line counts.
    :type max_sentinel: int
    :yield: A chunk of the file
    :rtype: Iterator[str]
    """
    cnt, chunk = 0, ''
    for line in f_obj:
        match = re.search(sentinel, line)
        if match:
            cnt += 1
        if cnt <= max_sentinel:
            chunk += line
        else:
            yield chunk
            cnt = 0
            chunk = line
    yield chunk

def process_chunk(this_chunk, my_sample_list):
    """Function to process a single line."""
    # my_dict: dict[str, list[tuple[str, int]]] = dict()
    my_dict = dict()

    lines = this_chunk.splitlines()
    for line in lines:
        if line.startswith('#'):
            pass
        else:
            # Split data lines into 10 fields, the last one is the samples info
            field_list = line.rstrip().split('\t', 9)
            chrom = field_list[0]
            pos = field_list[1]
            coord = '{}:{}'.format(chrom, pos)
            # Split the last field (sample info) and only keep the genotype
            for i, sample_info in enumerate(field_list[9].split('\t')):
                try:
                    # GT:DP:AD:GQ:GL
                    dp = int(sample_info.split(':')[1])
                except IndexError:
                    dp = 0
                if my_sample_list[i] not in my_dict:
                    my_dict[my_sample_list[i]] = dict()
                if chrom not in my_dict[my_sample_list[i]]:
                    my_dict[my_sample_list[i]][chrom] = dict()
                my_dict[my_sample_list[i]][chrom][pos] = dp
    return my_dict


# def parse_file_parallel(file_path, my_sample_list, num_workers=os.cpu_count()):
#     """Parses a file line by line in parallel."""
#     result_list = list()
#     with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
#         with open(file_path, 'r') as file:
#             # futures = [executor.submit(process_line, line, my_sample_list) for line in file]
#             # futures = [executor.submit(process_chunk, my_chunk, my_sample_list) for my_chunk in chunk_read(file, '\n', 1000)]
#             futures = [ executor.map(process_chunk, my_chunk, my_sample_list) for my_chunk in chunk_read(file, '\n', 1000)]
#             for future in concurrent.futures.as_completed(futures):
#                 result_list.append(future.result())
#
#     # Merge all results into a single dictionary
#     cov_dict: dict[str: int] = dict()
#     for result in result_list:
#         cov_dict = merge_dicts_append_lists([cov_dict, result])
#     return cov_dict

def parse_file_parallel(file_path, my_sample_list, num_workers=os.cpu_count()):
    cov_dict = dict()
    with futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        with open(file_path, 'r') as file:
            args = ((my_chunk, sample_list) for my_chunk in chunk_read(file, '\n', 1000))
            for result in executor.map(lambda x: process_chunk(*x), args):
                # cov_dict = merge_dicts_append_lists([cov_dict, result])
                cov_dict.update(result)
    return cov_dict

def merge_dicts_append_lists(dicts):
    merged = dict()
    for d in dicts:
        for key, value in d.items():
            if key in merged:
                if isinstance(merged[key], list):
                    if isinstance(value, list):
                        merged[key].extend(value)
                    else:
                        merged[key].append(value)
                else:
                    merged[key] = [merged[key], value]
            else:
                if isinstance(value, list):
                    merged[key] = value
                else:
                    merged[key] = [value]
    return merged

def coverage_graph(my_vcf_dict):
    # Coverage per sample
    df = pd.DataFrame(my_vcf_dict)
    df = df.astype(int)  # Convert type to integer
    # df = pd.DataFrame(columns=['Sample', 'DP'], index=None)
    # for sample, dp_list in my_vcf_dict.items():
    #     for i in range(len(dp_list)):
    #         dp = dp_list[i]
    #         data = {'Sample': [sample], 'DP': [dp]}  # Convert to dictionary
    #         tmp_df = pd.DataFrame.from_dict(data)  # Convert to dataframe
    #         df = pd.concat([df, tmp_df])  # Add at bottom of master dataframe

    # Remove missing and no coverage loci
    # df.replace(r'\.|\./\.', np.nan, regex=True, inplace=True)  # Replace missing by nan
    # df1 = df.loc[:, ['Sample', 'DP']].dropna()  # Extract depth and drop loci with no coverage
    # df1['DP'] = df1['DP'].astype(int)  # Convert type to integer
    # print(df1)

    #  Display a table summarizing the key statistical values of your boxplot.
    # stats_df = df1.groupby('Sample')['DP'].agg(
    #     Minimum='min',
    #     Q1=lambda x: x.quantile(0.25),
    #     Median='median',
    #     Mean='mean',
    #     Q3=lambda x: x.quantile(0.75),
    #     Maximum='max'
    # ).round(1).reset_index()

    # # Sort values ascending in dataframe
    # df.sort_index(ascending=True, inplace=True)

    # Split index column.
    df = df.reset_index(names='Coordinates')
    df[['chrom', 'pos']] = df['Coordinates'].str.split(':', expand=True, n=2)
    df['pos'] = df['pos'].astype(int)  # Convert column 'pos' to int
    df = df.iloc[:, 1:]  # Drop the first column (index 0)

    # Move the last two columns of a Pandas DataFrame to the first two positions,
    cols = df.columns.tolist()  # Get the list of all column names
    last_two_cols = cols[-2:]  # Extract the last two columns
    remaining_cols = cols[:-2]  # Get the remaining columns (all except the last two)
    new_col_order = last_two_cols + remaining_cols  # Create the new column order
    df = df[new_col_order]  # Reindex the DataFrame with the new column order

    df.sort_values(by=['chrom', 'pos'], ascending=[True, True], inplace=True)

    # Compute statistics
    stat_df = df.iloc[:, 2:]  # Drop the first two column (index 0, 1)
    summary_df = pd.DataFrame({
        'min': stat_df.min(),
        'Q1': stat_df.quantile(0.25),
        'Q2 (median)': stat_df.median(),
        'Q3': stat_df.quantile(0.75),
        'mean': stat_df.mean(),
        'max': stat_df.max()
    }).round(1)

    # Reorder index to be columns as rows
    summary_df = summary_df.T  # Optional: transpose if you want stats as rows
    #
    # title = 'Detailed of coverage statistics per sample: {}'.format(os.path.basename(my_vcf))
    #
    # fig_table = go.Figure(data=[go.Table(
    #     header=dict(values=list(df.columns)),
    #     cells=dict(values=[df[col].to_list() for col in df.columns])
    # )])
    # fig_table.update_layout(title_text=title)
    # fig_table.update_layout({'margin': {'t': 50}})

    # return fig_table, df, summary_df
    # return df, summary_df
    return df, summary_df

def make_plot(my_df, the_vcf):
    plt.figure(figsize=(16, 10))  # Width=12, Height=8 inches
    sns.boxplot(my_df, showfliers=False)
    # sns.violinplot(my_df, cut=0)
    plt.title('SNP Depth of Coverage Boxplot - No Outliers')
    plt.xlabel('Sample')
    plt.xticks(rotation=60)
    plt.ylabel('Depth of Coverage')

    out_plot = the_vcf.replace('.vcf', '_no_outliers.pdf')
    plt.savefig(out_plot, format='pdf', bbox_inches='tight')

    # Clear the current figure
    plt.clf()

    plt.figure(figsize=(16, 10))  # Width=12, Height=8 inches
    sns.boxplot(my_df)
    # sns.violinplot(my_df, cut=0)
    plt.title('SNP Depth of Coverage Boxplot - Outliers Included')
    plt.xlabel('Sample')
    plt.xticks(rotation=60)
    plt.ylabel('Depth of Coverage')

    out_plot = the_vcf.replace('.vcf', '_with_outliers.pdf')
    plt.savefig(out_plot, format='pdf', bbox_inches='tight')


# Run the stuff
if __name__ == "__main__":
    sample_list = get_sample_list(my_vcf)
    print('Parsing VCF file: {}'.format(os.path.basename(my_vcf)))
    coverage_dict = parse_vcf(my_vcf, sample_list)
    # coverage_dict = parse_file_parallel(my_vcf, sample_list)

    print('Generating coverage statistics...')
    # fig_table, df_stats, summary_df = coverage_graph(coverage_dict)
    df_stats, summary_df = coverage_graph(coverage_dict)

    # # Write to html
    # out_html = my_vcf.replace('.vcf', '.html')
    # with open(out_html, 'w') as f:
    #     f.write(fig_table.to_html(full_html=False, include_plotlyjs='cdn'))

    # Write to TSV
    out_table = my_vcf.replace('.vcf', '.tsv')
    pd.DataFrame.to_csv(df_stats, out_table, sep='\t',header=True, index=False)

    # Write to TSV
    out_table = my_vcf.replace('.vcf', '_summary.tsv')
    pd.DataFrame.to_csv(summary_df, out_table, sep='\t', header=True, index=True)

    make_plot(df_stats, my_vcf)
