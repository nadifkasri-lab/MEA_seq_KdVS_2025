#!/usr/bin/python3

"""
Script to store gene metadata information in a MongoDB database.

Example usage: python3 store_gene_info.py --file=GSE70138_Broad_LINCS_gene_info_2017-03-06.txt
					      --database=LINCS_phase2
				          --collection=gene_info

--database argument must be the authentication database of the user which is
assumed to have write rights to the db. If not, authentication upon logging in
to the MongoDB will fail.

If the collection already exists in the database, it will be replaced.
"""

import pandas as pd
from argparse import ArgumentParser
from pymongo import MongoClient
import sys
import numpy as np

parser = ArgumentParser()
parser.add_argument('-f', '--file', dest='gene_info_file',
        type=str, help='Path to gene information file', required=True)
parser.add_argument('-d', '--database', dest='destination_db',
        type=str, help='Database to store gene information in', required=True)
parser.add_argument('-c', '--collection', type=str, dest='destination_coll',
        help='Collection to store gene information in', required=True)


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

db = mongo_client[args.destination_db]

colls = db.list_collection_names()
print(f"Available collections in this database:{colls}")

if args.destination_coll in db.list_collection_names():
    print(
            f'Collection {args.destination_coll} already existed - dropping contents')
    getattr(db, args.destination_coll).remove()

else:
    print(f'Created new collection: {args.destination_coll}')

collection = db[args.destination_coll]

gene_info = pd.read_csv(args.gene_info_file, sep='\t')
gene_info.columns = ['gene_id', 'gene_symbol', 'gene_title', 'is_lm', 'is_bing']

gene_info['position'] = pd.Series(np.arange(1, gene_info.shape[0] + 1))

for row in gene_info.itertuples(index=False):
    row = row._asdict()

    res = collection.insert_one(row)

    if res.acknowledged == False:
        print(f'Insert failed: {row}')
        sys.exit(1)

print('All inserts succesfull')
