#!/usr/bin/python3

"""
Script to calculate the connectivity scores for user-specified up- and downregulated genes.

Genes must be supplied in .grp format:
More information: https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

User must have access to the MONGO database. Login.cfg file may be specified
on the command line. This config file must have a section LoginCredentials
which contains at least MongoLogin, MongoPassword, MongoHost and MongoPort:

[LoginCredentials]
MongoLogin=xxx
MongoPassword=xxx
MongoHost=cmbinas4.umcn.nl
MongoPort=27017

Gene expression profiles in the LINCS database are filtered in two separate
wyas: First, BING genes as defined by Cmap are removed. Furthermore, a file with
genes considered in the RNA-seq experiment may be specified using the --included
option. These are the genes that will be used to perform the calculations

The query log is generated in the output directory and contains information
about the amount of genes finally included.

Results are stored as a .pickle file: A compressed list of dictionaries, each
containing results for a single comparison to a gene expression profile.
Afterwards, run quantile.py to obtain the final results.

Example usage:

python3 query_controller.py
    --up up.grp
    --down down.grp
    --outdir ~/home/path/to/outdir
    --all
    --njobs 10
    --included genes.grp
    --login ~/login.cfg
"""

from query import Query, tau
import pandas as pd
import numpy as np
import pickle
from statistics import mean
from argparse import ArgumentParser
import sys
import os
import logging
import configparser


def parse_grp(filename: str) -> list:
    """
    Parse .grp file. Each line is assumed to contain a single entrez gene id
    Lines starting with # are ignored

    Arguments:
        filename (str): File to open

    Returns:
        List of integers containing entrez gene ids from file
    """

    genes = []

    try:
        with open(filename, 'r') as fh:
            lines = fh.readlines()

        for line in lines:
            if line.startswith('#'):
                continue
            else:
                genes.append(line.strip())

    except FileNotFoundError as e:
        logging.critical(f'Supplied .grp file not found: {filename}')
        sys.exit(1)

    except Exception as e:
        logging.critical(f'Error while parsing .grp file: {filename}: {e}')
        sys.exit(1)

    return [int(gene) for gene in genes]


if __name__ == "__main__":

    # Command line arguments
    parser = ArgumentParser(usage=__doc__)
    parser._action_groups.pop()

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--up', required=True, dest='up_genes',
    help='Path the .grp file containing upregulated genes')
    required.add_argument('--down', required=True, dest='down_genes',
    help='Path to the .grp file containing downregulated genes')
    required.add_argument('--included', required=True,
                        help = '.grp file that contains entrez gene ids of all the genes that were measured')
    required.add_argument('--login', required=True,
                        help = '(Default: login.cfg) Login configuration file. Is assumed to contain a \
                            single section LoginCredentials, \
                            containing MongoLogin, MongoPassword, MongoHost, MongoPort')


    # Optional
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--outdir', required=False, dest='outdir',
                        help='Output directory. Will be made if does not exist. Default: current directory')
    optional.add_argument('--all', required=False, dest='all', action='store_true',
                        help='If checked, searches all cellines and treatments. Default: True')
    optional.add_argument('--njobs', required=False, default=4, type=int,
                        help='Number of jobs using joblib backend for parallel calculations. Default: 4')
    args = parser.parse_args()

    # Check if the output directory exists. If not, make
    if args.outdir is not None:
        if not os.path.isdir(args.outdir):
            os.mkdir(args.outdir)

    # Logging. New log is created per query
    # If outdir argument is given, store log there. Else, in working dir.
    fname = os.path.join(
        args.outdir, 'query.log') if args.outdir is not None else 'query.log'

    # Configuration of the logging file
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-4s %(message)s',
                        datefmt='%m-%d %H:%M',
                        handlers=[logging.FileHandler(filename=fname, mode='w')])

    logging.info(f'Reading upregulated genes from: {args.up_genes}')
    logging.info(f'Reading downregulated genes from: {args.down_genes}')

    if not args.up_genes.endswith('.grp') or not args.down_genes.endswith('.grp'):
        logging.critical('Supplied genes must be in .grp format')
        sys.exit(1)

    # Open query genes and included genes.
    up_genes = parse_grp(args.up_genes)
    down_genes = parse_grp(args.down_genes)
    included = parse_grp(args.included)

    logging.info(f'Total up genes: {len(up_genes)}')
    logging.info(f'Total down genes: {len(down_genes)}')
    logging.info(f'Total included genes: {len(included)}')

    # Cellines and treatments to search in
    # If unspecified, we search only our own cellines
    if args.all == False:
        cellines_phase1 = ['FIBRNPC', 'NPC', 'VCAP',
                           'A375', 'A549', 'SHSY5Y', 'NEU', 'HEK293T']
        cellines_phase2 = ['NPC', 'NPC.CAS9', 'NPC.TAK',
                           'A375', 'A549', 'NEU', 'MNEU.E', 'HUES3']

        treatment_phase1 = ['trt_cp', 'trt_sh.css', 'trt_sh.cgs', 'trt_sh']
        treatment_phase2 = ['trt_cp']

    # Else we search all cellines in the database
    else:
        phase_1 = pd.read_csv(
            '~/LINCS_data/phase_1/metadata/GSE92742_Broad_LINCS_sig_info.txt', sep='\t', low_memory=False)
        cellines_phase1 = phase_1.cell_id.unique()
        treatment_phase1 = phase_1.pert_type.unique()


        phase_2 = pd.read_csv(
            '~/LINCS_data/phase_2/metadata/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt', sep='\t')
        cellines_phase2 = phase_2.cell_id.unique()
        treatment_phase2 = phase_2.pert_type.unique()

    logging.info('Searching cellines:')
    logging.info(f"Phase 1: {' - '.join(cellines_phase1)}")
    logging.info(f"Phase 2: {' - '.join(cellines_phase2)}")
    logging.info('Searching treatments: ')
    logging.info(f"Phase 1: {' - '.join(treatment_phase1)}")
    logging.info(f"Phase 2: {' - '.join(treatment_phase2)}")

    # Import the Mongo credentials
    try:
        config = configparser.ConfigParser()

        logging.info(f'Reading config file: {args.login}')
        config.read(args.login)

        mongo_credentials = {
            'USERNAME': config.get('LoginCredentials', 'MongoLogin'),
            'PASSWORD': config.get('LoginCredentials', 'MongoPassword'),
            'MONGO_HOST': config.get('LoginCredentials', 'MongoHost'),
            'MONGO_PORT': config.getint('LoginCredentials', 'MongoPort')
        }
    except Exception as e:
        logging.critical(f'Error while reading login file: {e}')

    # Actual querying database and calculating scores
    query_phase1 = Query(database='LINCS_phase1',
                         collection='data',
                         mongo_credentials=mongo_credentials,
                         njobs=args.njobs)

    logging.info('Starting search in LINCS_phase1 database')

    results_phase1 = query_phase1.search_db(cellines=cellines_phase1,
                                            treatments=treatment_phase1,
                                            query_up=up_genes,
                                            query_down=down_genes,
                                            included = included)

    results_phase1 = query_phase1.normalize_wcs(results_phase1)

    logging.info(
        f'Done with search in LINCS_phase1 database. {len(results_phase1)} results')
    logging.info('Starting search in LINCS_phase2 database')

    query_phase2 = Query(database='LINCS_phase2',
                         collection='data',
                         mongo_credentials=mongo_credentials,
                         njobs=args.njobs)

    results_phase2 = query_phase2.search_db(cellines=cellines_phase2,
                                            treatments=treatment_phase2,
                                            query_up=up_genes,
                                            query_down=down_genes,
                                            included = included)

    results_phase2 = query_phase2.normalize_wcs(results_phase2)

    logging.info(
        f'Done with search in LINCS_phase2 database. {len(results_phase2)} results')

    # Combine results and calculate tau score
    results = results_phase1 + results_phase2
    results = tau(results)
    logging.info(f'Pickling: {len(results)} results')

    # Save results in specified output dir if supplied
    fname = os.path.join(
        args.outdir, 'results.pickle') if args.outdir is not None else 'results.pickle'

    logging.info(f'Storing results in {fname}')

    with open(fname, 'wb') as fh:
        pickle.dump(results, fh)
