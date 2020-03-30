import json
from pathlib import Path

import pandas as pd


def read_input(quant_dirs):
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


def main(quant_dirs, accession_map, fname_out):
    # read quantification results
    df = read_input(quant_dirs)

    # add additional information
    df['name'] = df['accession'].map(accession_map)

    df.set_index(['name', 'target_id'], inplace=True)

    df['relative_count'] = df['est_counts'] / df['read_number']

    df['relative_count_fraction'] = df.groupby('target_id')['relative_count'].transform(lambda x: x / x.loc['AX_2_WT'].iloc[0])

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
        snakemake.input.quant_dirs,
        snakemake.config['samples'],
        snakemake.output.fname)
