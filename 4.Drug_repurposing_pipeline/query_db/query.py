#!/usr/bin/python3

"""
Class to perform a query in MongoDB servers with user-supplied up and down-regulated genes

Do not call this script directly, rather use it through query_controller.py
"""


from pymongo import MongoClient
from joblib import Parallel, delayed
import itertools as it
import sys
import numpy as np
import pandas as pd
from statistics import mean
import logging

class Query:
    """
    Class to search MongoDB, perform calculations in parallel and aggregate the results.
    Based on supplied cell lines and treatment options in the, a joblib backend is started with
    default 4 jobs to allow parallel retrieval of data from database and calculations.

    Call and use this class through the controller script.
    """

    def __init__(self, database: str, collection: str, mongo_credentials: dict, njobs: int):

        self.database = database
        self.collection = collection
        self.mongo_credentials = mongo_credentials
        self.nr_jobs = njobs

        # Test if mongoDB credentials are valid
        _ = self.login_mongo(credentials=self.mongo_credentials)
        logging.info('Mongo connection succesfull')


    def login_mongo(self, credentials: dict) -> MongoClient:
        """Helper function to login to Mongo and return the client if connection is succesfull

        Args:
          credentials (dict): Contains at least the fields MONGO_HOST, MONGO_PORT, USERNAME, PASSWORD.
          Other fields will be ignored.

        Returns:
          pymongo::MongoClient instance
        """

        try:
            client = MongoClient(host = credentials['MONGO_HOST'],
                                 username = credentials['USERNAME'],
                                 password = credentials['PASSWORD'],
                                 port = credentials['MONGO_PORT'],
                                 authSource = self.database)

            # Will throw error if connection cannot be made
            _ = client.server_info();

            return client

        except KeyError as e:
            logging.critical('Supplied credentials do not contain required fields')
            sys.exit(1)

        except Exception as e:
            logging.critical(f'An error has ocurred while logging in to Mongo: {e}')
            sys.exit(1)

    def search(self, cell_id: str, treatment: str) -> list:
        """
        Single search in the MongoDB based on specific cell_id / treatment

        Args:
            cell_id: particular cell_id to perform search on
            treatments: partcular treatment to perform search on

        Returns:
            List of dictionaries. Each dict contains a calculation against a
            database signature. Dictionaries contain fields 'metadata', and
            'wcs'
        """
        client = self.login_mongo(credentials=self.mongo_credentials)
        coll = client[self.database][self.collection]

        filter = {'metadata.cell_id': cell_id,
                'metadata.pert_type': treatment}

        # Search MongoDB but excluse _id field from return values
        query_res = coll.find(filter=filter, projection={'_id': False})

        # Calculate the connectivity score for this portion of the results
        conn_scores = self.weighted_conn_score(query_res)

        return conn_scores


    def search_db(self, cellines: list, treatments: list, query_up: list, query_down: list, included: list) -> list:
        """
        Search MongoDB based on the supplied cellines and treatments. Every combination between these two lists
        is produced and searched for. Log files contain the number of documents returned per search.

        Args:
            cellines: List of cellines (strings) to be searched for in MongoDB
            Treatments: List of treatments (strings) to be searched for in MongoDB
            query_up: List containing sorted upregulated genes - First gene in the array is most
                differentially expressed
            queyr_down: List containing downregulated genes - First gene in the array is most
                differentially expressed

        Returns:
            List of dictionaries - 1 dict per search performed.
        """

        # Find how many of the genes are measured by LINCS and thus valid
        client = self.login_mongo(credentials=self.mongo_credentials)
        coll = client[self.database]['gene_info']

        cursor = coll.find(filter = {}, projection = {'_id': False})
        df = pd.DataFrame.from_dict({pos: row for pos, row in enumerate(cursor)}, 'index')
        
        # We store the overlap between the genes that we include and the genes measured by LINCS
        self.included_genes = [gene for gene in included if gene in df.gene_id.to_list()]
        logging.info(f'Overlap between included genes and LINCS genes: {len(self.included_genes)}')
        
        # Subset only the BING gene space to use for the calculations
        logging.info('Validating input genes. Using only BING genes')
        bing_genes = df[df['is_bing'] == 1]['gene_id'].to_list()

        # Filter up and down genes for the proper genes that we know are included
        query_up_filtered = [item for item in query_up if item in bing_genes]
        query_down_filtered = [item for item in query_down if item in bing_genes]

        logging.info(f'Retaining: {len(query_up_filtered)} up-regulated genes out of {len(query_up)}')
        logging.info(f'Retaining: {len(query_down_filtered)} down-regulated genes out of {len(query_down)}')

        self.query_up = query_up_filtered
        self.query_down = query_down_filtered

        # All combinations of the cell lines and treatments
        combs = it.product(cellines, treatments)
        coll = client[self.database]['data']

        total = 0
        for cell_id, trt_type in it.product(cellines, treatments):
            count = coll.find({'metadata.cell_id': cell_id, 'metadata.pert_type': trt_type}).count()
            # logging.info(f'Celline: {cell_id} - Treatment: {trt_type}: {count} perturbation profiles in db')
            total += count

        logging.info(f'Total results in database {self.database}: {total}')

        # Search parallel using joblib
        logging.info(f'Searching database using {self.nr_jobs} parallel threads')

        total = Parallel(n_jobs = self.nr_jobs)(delayed(self.search)(*params) for params in combs)

        # Total is a list of lists that contains the results per worker. Flatten list and return
        return [item for sublist in total for item in sublist]


    def weighted_conn_score(self, cursor) -> list:
        """
        Calculates the weighted connectivity score as explained in Subramanian et al (2017)

        Args:
            cursor returned by the MongoDB search. Contains results for a specific cellid / treatment combination

        Returns:
            List of dictionaries. 1 dict per search performed.
        """

        def KS(query: list, df: pd.DataFrame) -> float:
            """
            Helper function to calculate KS statistic

            Args:
                query: np.array containing either up or down-regulated query genes
                df: contains gene_symbols column - 12328 gene symbols in right order
            """

            # Sort df on decreasing z-scores
            df = df.sort_values(by=['z_scores'], ascending=False)
            df.reset_index(inplace=True, drop=True)
            
            # n is the length of the gene space we search in
            n = df.shape[0]

            # Subset df to only retain query genes
            df = df[df['gene_id'].isin(query)]
            indices = df.index.to_numpy()

            # t is the length of the query 
            t = len(indices)
            
            j = np.arange(1, t + 1)
            a = np.max((j / t) - (indices / n))

            j = np.arange(0, t)
            b = np.max((indices / n) - (j / t))

            return (a if a > b else -b)


        # Pull gene information from the database
        client = self.login_mongo(credentials = self.mongo_credentials)
        coll = client[self.database]['gene_info']

        gene_meta = coll.find({}, {'_id': False})
        df = pd.DataFrame.from_dict({pos: row for pos, row in enumerate(gene_meta)}, 'index')

        total = []
        for doc in cursor:

            result = {}

            result['metadata'] = doc['metadata']
            df['z_scores'] = doc['data']
            
            # Take out the non-BING genes
            bing = df[df['is_bing'] == 1]
          
            # Take only the genes that we have measured as well  
            final = bing[bing['gene_id'].isin(self.included_genes)]
            
            ks_up = KS(self.query_up, final)
            ks_down = KS(self.query_down, final)

            if (ks_up < 0) and (ks_down < 0):
                result['wcs'] = 0.0
            elif (ks_up > 0) and (ks_down > 0):
                result['wcs'] = 0.0
            else:
                result['wcs'] = (ks_up - ks_down)

            total.append(result)

        return total


    def normalize_wcs(self, wcs_list: list) -> list:
        """ Normalizes the wcs scores in given list
        according to Subramanian et al (2017)

        Args:
            wcs_list: List of dictionaries that contains the following fields:
                metadata and wcs.

        Returns:
            List of dictionaries identical to wcs_list, with added 'ncs' entry
            per dict
        """

        # Get all unique cell lines and perturbants
        cell_lines = list(set(d['metadata']['cell_id'] for d in wcs_list))
        perturbants = list(set(d['metadata']['pert_type'] for d in wcs_list))

        result = []

        for cell_line in cell_lines:

            for perturbant in perturbants:

                # Select only those results in the list that have this cell_line and
                # pert type, separate + and -

                pos_mu = [item for item in wcs_list if
                          (item['metadata']['cell_id'] == cell_line) and
                          (item['metadata']['pert_type'] == perturbant) and
                          (item['wcs'] > 0.0)]

                neg_mu = [item for item in wcs_list if
                          (item['metadata']['cell_id'] == cell_line) and
                          (item['metadata']['pert_type'] == perturbant) and
                          (item['wcs'] < 0.0)]

                # Only calculate the means if there is acutally something in the list
                if pos_mu:
                    mean_pos = mean([item['wcs'] for item in pos_mu])

                    for item in pos_mu:
                        item['ncs'] = {}
                        item['ncs']['score'] = item['wcs'] / mean_pos
                        item['ncs']['sign'] = 'positive'

                    result.append(pos_mu)


                if neg_mu:
                    mean_neg = mean([item['wcs'] for item in neg_mu])

                    for item in neg_mu:
                        item['ncs'] = {}
                        item['ncs']['score'] = item['wcs'] / mean_neg
                        item['ncs']['sign'] = 'negative'

                    result.append(neg_mu)

                else:
                    pass
                    # print(f'No results for: {cell_line} - {perturbant}')

        # Flatten list of results and return
        return [item for sublist in result for item in sublist]


def tau(ncs_list: list) -> list:
    """
    Calculates tau score as described in Subramanian et al 2017

    Args:
        ncs_list: List of dictionaries

    Returns:
        List of dictionaries identical to ncs_list, with added 'tau' field
    """

    scores = []
    n = len(ncs_list)

    for item in ncs_list:
        if item['ncs']['sign'] == 'positive':
            scores.append(item['ncs']['score'])
        else:
            scores.append(-1 * item['ncs']['score'])

    scores = np.array(scores)
    scores_abs = np.absolute(scores)

    for pos, item in enumerate(ncs_list):
        smaller = (scores_abs < abs(item['ncs']['score'])).sum()

        tau = (100 / n) * smaller

        if item['ncs']['sign'] == 'negative':
            tau = -1 * tau

        item['tau'] = tau

        ncs_list[pos] = item

    return ncs_list


if __name__ == "__main__":

    pass
