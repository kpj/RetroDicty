import json
from pathlib import Path

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import pysam
from Bio import SeqIO


sns.set_context('talk')


COVERAGEPLOT_YLIMITS = {
    'DIRS-1': (-35_000, 20_000),
    'Skipper-1_PLASMID': (-4_500, 2_000)
}


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


def coverage_plot(accession, ref, coverage, coverage_rev, mapping_count, fname):
    plt.figure(figsize=(8, 6))

    sns.lineplot(x=np.arange(0, len(coverage)), y=coverage, label='forward')
    sns.lineplot(x=np.arange(0, len(coverage_rev)), y=-coverage_rev.astype(int), label='reverse')

    plt.fill_between(np.arange(0, len(coverage)), coverage)
    plt.fill_between(np.arange(0, len(coverage_rev)), -coverage_rev.astype(int))

    plt.title(f'{accession} - {ref} (Mapped reads: {mapping_count})')
    plt.xlabel('Position [bp]')
    plt.ylabel('Read coverage')

    plt.legend(loc='best')

    if ref in COVERAGEPLOT_YLIMITS:
        plt.ylim(COVERAGEPLOT_YLIMITS[ref][0], COVERAGEPLOT_YLIMITS[ref][1])

    plt.tight_layout()
    plt.savefig(fname)


def count_fastq_reads(fname):
    with open(fname) as fd:
        return int(sum(1 for line in fd) / 4)


def get_read_callback(reverse=False):
    def tmp(read):
        if read.is_unmapped or read.is_secondary or read.is_qcfail or read.is_duplicate:
            return False

        return read.is_reverse == reverse
    return tmp


def read_input_star(fastq_files, bam_files, fname_ref, plot_dir):
    tmp = []

    # parse reference sequences
    ref_seqs = {}
    for record in SeqIO.parse(fname_ref, 'fasta'):
        ref_seqs[record.id] = str(record.seq)

    # parse mapping results
    for fname_fastq, fname_bam in zip(fastq_files, bam_files):
        bam = pysam.AlignmentFile(fname_bam)
        assert bam.has_index()

        accession = str(fname_bam).split('/')[1]
        read_number = count_fastq_reads(fname_fastq)

        for ref in bam.references:
            mapping_count = bam.count(contig=ref)

            coverage = np.sum(
                bam.count_coverage(
                    contig=ref,
                    read_callback=get_read_callback(reverse=False))
                , axis=0)
            coverage_rev = np.sum(
                bam.count_coverage(
                    contig=ref,
                    read_callback=get_read_callback(reverse=True))
                , axis=0)

            # coverage plot
            coverage_plot(
                accession, ref,
                coverage, coverage_rev,
                mapping_count,
                plot_dir / f'coverage_{accession}_{ref}.pdf')

            # store coverage data
            with open(plot_dir / f'coverage_{accession}_{ref}.csv', 'w') as fd:
                fd.write(
                    'position,' +
                    ','.join(str(p) for p in range(len(ref_seqs[ref]))) +
                    '\n')
                fd.write('reference,' + ','.join(ref_seqs[ref]) + '\n')
                fd.write(
                    'coverage_forward,' +
                    ','.join(str(c) for c in coverage) +
                    '\n')
                fd.write(
                    'coverage_reverse,' +
                    ','.join(str(c) for c in coverage_rev) +
                    '\n')

            # store data
            tmp.append({
                'accession': accession,
                'reference': ref,
                'total_read_number': read_number,
                'mapped_read_number': mapping_count
            })

    return pd.DataFrame(tmp)


def main(fastq_files, bam_files, fname_ref, accession_map, fname_out, plot_dir):
    # read quantification results
    plot_dir.mkdir(exist_ok=True)
    df = read_input_star(fastq_files, bam_files, fname_ref, plot_dir)

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
        snakemake.input.fname_ref,
        snakemake.config['samples'],
        snakemake.output.fname, Path(snakemake.output.plot_dir))
