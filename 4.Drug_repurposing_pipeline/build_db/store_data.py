#!/usr/bin/python3

"""
Store perturbation experiment results in MongoDB (lists of z-scores).

Example usage: python3 store_data.py --file GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
                                     --metadata GSE92742_Broad_LINCS_sig_info.txt
                                     --database LINCS_phase1
                                     --collection data

--database argument must be the authentication database else login to MongoDB
will fail

If a collection with the same name already exists it will be dropped.
"""

import pandas as pd
from argparse import ArgumentParser
from pymongo import MongoClient
import sys
import numpy as np
from cmapPy.pandasGEXpress.parse import parse

parser = ArgumentParser()

parser.add_argument('-f', '--file', dest='gctx_file',
                    type=str, help='Path to .gctx file', required=True)
parser.add_argument('-m', '--metadata', dest='sig_info',
                    help='Sig_info metadata file', type=str, required=True)
parser.add_argument('-d', '--database', dest='destination_db', type=str,
                    help='Database to store information in', required=True)
parser.add_argument('-c', '--collection', dest='destination_coll', type=str,
                    help='Collection to store information in', default='data')

args = parser.parse_args()


MONGO_HOST = 'cmbinas4.umcn.nl'
MONGO_PORT = 27017
USERNAME = 'martijnz'
PASSWORD = 'Z990141'

try:
    mongo_client = MongoClient(host=MONGO_HOST,
            username=USERNAME,
            password=PASSWORD,
            authSource=args.destination_db,
            port=MONGO_PORT)

    # If the client is not connected, this will throw an error
    mongo_client.server_info()

except Exception as e:
    print(e)
    sys.exit(1)

# Create the list from which to sort the df
coll = mongo_client['LINCS_phase1']['gene_info']
cursor = coll.find(filter = {}, projection = {'_id': False})
meta = pd.DataFrame.from_dict({pos: row for pos, row in enumerate(cursor)}, 'index')

meta.set_index(meta['gene_id'], inplace=True, drop=True)
new_idx = meta.index

db = mongo_client[args.destination_db]

print(
    f'Storing in database: {args.destination_db} - collection: {args.destination_coll}')

if args.destination_coll in db.list_collection_names():
    print(
        f'Collection {args.destination_coll} already existed - dropping contents')
    getattr(db, args.destination_coll).remove()

else:
    print(f'Created new collection: {args.destination_coll}')

collection = db[args.destination_coll]

# Parse each celline individually
sig_info = pd.read_csv(args.sig_info, sep='\t')
cellines = sig_info["cell_id"].unique()

total = len(cellines)

for pos, celline in enumerate(cellines):
    print(f'Celline: {pos} out of {total}')

    coldata = sig_info[sig_info["cell_id"] == celline]

    identifiers = coldata["sig_id"]

    # Parsing of actual .gctx file - only specified celline
    gctx = parse(args.gctx_file, cid=identifiers)

    for row in coldata.itertuples(index=False):

        # Final object to be stored in the db
        complete = {}

        # Metadata about the experiment
        metadata = row._asdict()

        # LINCS_phase1 contains non-ascii characters - mu. Here, we check
        # for those and replace them with simple u.
        if all(ord(char) < 128 for char in metadata['pert_dose_unit']):
            pass

        # The only non-ascii char is mu, we replace this with u.
        else:
            st = metadata['pert_dose_unit']
            st = ''.join([i if ord(i) < 128 else 'u' for i in st])
            metadata['pert_dose_unit'] = st

            st = metadata['pert_idose']
            st = ''.join([i if ord(i) < 128 else 'u' for i in st])
            metadata['pert_idose'] = st

        # Actual z-scores
        z_scores = gctx.data_df[metadata['sig_id']]
        z_scores.index = z_scores.index.astype(int)

        before = z_scores.loc[5720]

        z_scores = z_scores.reindex(new_idx)
        after = z_scores.loc[5720]
        assert before == after

        z_scores = z_scores.to_list()

        complete['metadata'] = metadata
        complete['data'] = z_scores

        res = collection.insert_one(complete)

        if res.acknowledged:
            pass
        else:
            print('Query has not been acknowledged')
            sys.exit(1)
