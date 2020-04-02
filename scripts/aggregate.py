import json
from pathlib import Path

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import pysam
from Bio import SeqIO


sns.set_context('talk')


def read_input_kallisto(quant_dirs):
    df_list = []

    for dname in quant_dirs:
        dname = Path(dname)

        fname_abd = dname / 'abundance.tsv'
        tmp = pd.read_csv(fname_abd, sep='\t')

        fname_meta = dname / 'run_info.json'
        with open(fname_meta) as fd:
            meta = json.load(fd)
        tmp['read_number'] = meta['n_processed']

        accession = str(dname).split('/')[1]
        tmp['accession'] = accession

        df_list.append(tmp)

    return pd.concat(df_list)


def coverage_plot(accession, ref, coverage, mapping_count, fname):
    plt.figure(figsize=(8, 6))

    sns.lineplot(x=np.arange(0, len(coverage)), y=coverage)

    plt.title(f'{accession} - {ref} (Mapped reads: {mapping_count})')
    plt.xlabel('Position [bp]')
    plt.ylabel('Read coverage')

    plt.tight_layout()
    plt.savefig(fname)


def read_input_star(fastq_files, bam_files, plot_dir):
    tmp = []

    for fname_fastq, fname_bam in zip(fastq_files, bam_files):
        bam = pysam.AlignmentFile(fname_bam)
        assert bam.has_index()

        accession = str(fname_bam).split('/')[1]
        read_number = 3#len(list(SeqIO.parse(fname_fastq, 'fastq')))

        for ref in bam.references:
            mapping_count = bam.count(contig=ref)

            coverage_base = bam.count_coverage(contig=ref)
            coverage = np.sum(coverage_base, axis=0)

            # coverage plot
            coverage_plot(
                accession, ref,
                coverage, mapping_count,
                plot_dir / f'coverage_{accession}_{ref}.pdf')

            # store data
            tmp.append({
                'accession': accession,
                'reference': ref,
                'total_read_number': read_number,
                'mapped_read_number': mapping_count
            })

    return pd.DataFrame(tmp)


def main(fastq_files, bam_files, accession_map, fname_out, plot_dir):
    # read quantification results
    plot_dir.mkdir(exist_ok=True)
    df = read_input_star(fastq_files, bam_files, plot_dir)

    # add additional information
    df['name'] = df['accession'].map(accession_map)

    df.set_index(['name', 'reference'], inplace=True)

    df['relative_count'] = df['mapped_read_number'] / df['total_read_number']

    df['relative_count_fraction'] = df.groupby('reference')['relative_count'].transform(lambda x: x / x.loc['AX_2_WT'].iloc[0])

    df.reset_index(inplace=True)

    # pivot results
    # df_table = df.pivot(
    #     index='target_id',
    #     columns='name',
    #     values=['read_number', 'est_counts', 'relative_count', 'relative_count_fraction'])

    # save results
    df.to_csv(fname_out, index=False)


if __name__ == '__main__':
    main(
        snakemake.input.fastq_files, snakemake.input.bam_files,
        snakemake.config['samples'],
        snakemake.output.fname, Path(snakemake.output.plot_dir))
