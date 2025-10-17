#!/usr/bin/python3

"""
Calculates a cell-summarised score per compound.
Algorithm from: https://clue.io/connectopedia/cmap_algorithms

Example usage: python3 quantiles.py --pickle results.pickle
                                    --outdir results/
                                    --min_n_sig 10
"""

import pickle
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import os

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pickle', required=True, help='Input pickle as produced by query scripts')
    parser.add_argument('--outdir', required=True, help='Existing dir to store output')
    parser.add_argument('--min_n_sig', required=False, default=10, type=int, help = '(Default: 10). Minimum number of signatures to filter for per perturbagen')
    parser.add_argument('--filter', action='store_true', help = '(Default: False). Filter for only trt_cp instances')

    args = parser.parse_args()

    # Using same values as Cmap
    qhi = 67
    qlow = 33

    with open(args.pickle, 'rb') as fh:
        results = pickle.load(fh)

    # Filter out only compound treatments if specified
    if args.filter:
        results = [item for item in results if item['metadata']['pert_type'] == 'trt_cp']

    drugs = [item['metadata']['pert_iname'] for item in results]
    unique_drugs = list(set(drugs))

    # Convert the list of results to a dict for faster lookup
    # res is a dict that contains as keys the pert_iname, and as values
    # a list of dictionaries. Each dict in the value list contains the
    # contents of the metadata dict for this perturbation, along with the NCS
    res = {}
    for item in results:
        key = item['metadata']['pert_iname']
        if key not in res:
            res[key] = []

        res[key].append(item['metadata'])

        if item['ncs']['sign'] == 'positive':
            res[key][-1]['ncs'] = item['ncs']['score']
        else:
            res[key][-1]['ncs'] = -1 * item['ncs']['score']



    # Quantiles is the dict that will contain the final results for the quantile
    # calculations. In the case of multiple treatments for a single pert_iname
    # e.g. gene oe and trt_sh - keys are pert_iname _ pert_type. In the case
    # of a single treatment, key is pert_iname
    # Values are dicts that contain the quantile score, number of signatures and
    # the compound name for merging with annotation df
    quantiles = {}
    for drug in unique_drugs:

        # Subset contains as key the pert_iname, value is a list of
        # metadata dictionaries for each of the signatures
        subset = {k: v for k, v in res.items() if k == drug}

        # Get list of unique treatments for this perturbagen
        treatments = list(set([item['pert_type'] for item in list(*subset.values())]))

        for treatment in treatments:
            if len(treatments) == 1:
                name = drug
            else:
                name = f'{drug}_{treatment}'

            # Take the subset of the perturbation profiles for this specific pert_type
            sub = {}
            for key, value in subset.items():
                for item in value:
                    if item['pert_type'] == treatment:

                        if key not in sub:
                            sub[key] = []

                        sub[key].append(item)

            # Get the amount of signatures in this subset
            signatures = list(sub.values())[0]
            n = len(signatures)

            # Filter out perturbagens that have less signatures than the threshold
            if n < args.min_n_sig:
                continue

            n_pos = len([item for item in signatures if item['ncs'] > 0])
            n_neg = len([item for item in signatures if item['ncs'] < 0])

            scores = []
            for k, v in sub.items():
                for x in v:
                    scores.append(x['ncs'])

            scores = np.array(scores)

            # if the number of signatures is not enough, skip perturbagen

            upper = np.quantile(scores, qhi / 100)
            lower = np.quantile(scores, qlow / 100)

            if abs(upper) >= abs(lower):
                quantile = upper
            else:
                quantile = lower

            quantiles[name] = {}
            quantiles[name]['quantile'] = round(quantile, 3)
            quantiles[name]['quant_high'] = round(upper, 3)
            quantiles[name]['quant_low'] = round(lower, 3)
            quantiles[name]['n_sig'] = len(scores)
            quantiles[name]['n_sig_pos'] = n_pos
            quantiles[name]['n_sig_neg'] = n_neg
            quantiles[name]['trt'] = treatment
            quantiles[name]['pert_iname'] = drug

    df = pd.DataFrame.from_dict(quantiles, 'index')

    # Merge together with the annotation df
    annot = pd.read_csv('/home/martijnz/LINCS_data/drug_annotation/drug_annotation.csv')

    annot.set_index('drug', inplace=True, drop=True)
    annot = annot[annot.index.isin(df.index)]

    df = pd.merge(df, annot, how='outer', left_on = 'pert_iname', right_index = True)
    df.sort_values('quantile', inplace=True, ascending=False)
    df = df[['pert_iname', 'trt', 'quantile', 'quant_low', 'quant_high', 'n_sig', 'n_sig_pos', 'n_sig_neg', 'moa']]

    df.reset_index(inplace=True, drop=True)

    cols = ['sig_id', 'pert_id', 'pert_iname', 'cell_id', 'pert_type']
    phase_1 = pd.read_csv(
        '/home/martijnz/LINCS_data/phase_1/metadata/GSE92742_Broad_LINCS_sig_info.txt', sep='\t', low_memory=False)
    phase_2 = pd.read_csv(
        '/home/martijnz/LINCS_data/phase_2/metadata/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt', sep='\t', low_memory=False)
    meta = pd.concat([phase_1[cols], phase_2[cols]], axis=0)

    sizes = meta.groupby(['pert_iname', 'pert_type']).size().reset_index()
    sizes.columns = ['pert_iname', 'pert_type', 'n_sig_initial']

    # left join with the df - preserve key order in df
    df = pd.merge(df, sizes, how = 'left', left_on = ['pert_iname', 'trt'], right_on = ['pert_iname', 'pert_type'])
    df = df[['pert_iname', 'trt', 'quantile', 'quant_low', 'quant_high', 'n_sig', 'n_sig_pos', 'n_sig_neg', 'n_sig_initial', 'moa']]

    df.to_csv(os.path.join(args.outdir, 'quantiles.csv'))
