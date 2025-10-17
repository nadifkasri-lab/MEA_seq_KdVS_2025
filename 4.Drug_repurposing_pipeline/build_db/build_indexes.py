#!/usr/bin/python3


"""
Build indexes in the MongoDB.

Drops all existing indexes in the database

"""

import pandas as pd
from pymongo import MongoClient, ASCENDING, DESCENDING
import sys
import numpy as np
from cmapPy.pandasGEXpress.parse import parse

MONGO_HOST = 'cmbinas4.umcn.nl'
MONGO_PORT = 27017
USERNAME = 'martijnz'
PASSWORD = 'Z990141'

DATABASE = 'LINCS_phase1'

try:
    mongo_client = MongoClient(host=MONGO_HOST,
                        username=USERNAME,
                        password=PASSWORD,
                        authSource=DATABASE,
                        port=MONGO_PORT)

        # If the client is not connected, this will throw an error
    mongo_client.server_info()

except Exception as e:
        print(e)

print(f'Connected to database: {DATABASE}')

db = mongo_client[DATABASE]
coll = db['data']

coll.drop_indexes()

_ = coll.create_index([('metadata.cell_id', ASCENDING)], name="cell_id_index")
_ = coll.create_index([('metadata.pert_type', ASCENDING)], name='pert_type_index')
_ = coll.create_index([('metadata.sig_id', ASCENDING)], name='sig_id_index')
_ = coll.create_index([
                ('metadata.cell_id', ASCENDING),
                ('metadata.pert_type', ASCENDING),
                ], name='cell_trttype_dualindex')

DATABASE='LINCS_phase2'

try:
    mongo_client = MongoClient(host=MONGO_HOST,
                        username=USERNAME,
                        password=PASSWORD,
                        authSource=DATABASE,
                        port=MONGO_PORT)

        # If the client is not connected, this will throw an error
    mongo_client.server_info()

except Exception as e:
        print(e)
print(f'Connected to database: {DATABASE}')

db = mongo_client[DATABASE]
coll = db['data']

# Drop all indexes first - except the standard on _id field
coll.drop_indexes()


_ = coll.create_index([('metadata.cell_id', ASCENDING)], name='cell_id_index')
_ = coll.create_index([('metadata.pert_type', ASCENDING)], name='pert_type_index')
_ = coll.create_index([('metadata.sig_id', ASCENDING)], name='sig_id_index')
_ = coll.create_index([
                ('metadata.cell_id', ASCENDING),
                ('metadata.pert_type', ASCENDING),
                ], name='cell_trttype_dualindex')
