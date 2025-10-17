
These scripts implement the algorithms as described by Subramanian et al (2017). In short, users supply a gene expression profile in the form of up- and downregulated genes in .grp format (see: https://clue.io/connectopedia/grp_gmt_gmx_format).

Drug-perturbation profiles from the LINCS consortium are stored on internal MongoDB servers. This pipeline will compute the similarity / dissimilarity between the query gene expression profile and each gene expression profile in the database. More information about the algorithms: Subramanian et al 2017, Lamb et al 2006.


# Usage

Install required packages:

``` 
python3 -m pip install -r requirements.txt
``` 

An example script to run the relevant python scripts is included in the main directory: run_query.sh.
This script calls two python3 scripts to calculate the connectivity score per comparison between query and gene expression profile, and afterwards to perform
quantile normalisation.

Users can change the directory where all their data is stored directly in the shell script and submit to the cluster as:

``` 
qsub run_query.sh
``` 

# Output

The scripts produce two separate outputs: a compressed pickle file that contains the metadata and individual scores per signature in the reference database, and a comma-separated file that contains the final quantile-normalised scores per perturbagen.

The compressed pickle file may be opened with python as:

``` 
#!/bin/python3

import pickle

with open('results.pickle', 'rb') as fh:
    results = pickle.load(fh)

``` 
We can sort this list of dictionaries based on their associated tau score in increasing order as follows:

```
results = sorted(results, key = lambda x:x['tau'], reverse=False)

```


The pickle file contains a list of dictionaries. Each dictionary in the list belongs to a single comparison made between the query genes and the LINCS instances in the reference database. We removed comparisons that did not show significant connectivity to the query genes (e.g. WCS == 0.0), therefore the number of comparisons will differ per query.
A dictionary in the list contains the following information, where WCS is the weighted connectivity score and NCS is the normalised connectivity score.

```
{'metadata': {'cell_id': 'HA1E',
              'distil_id': 'REP.A025_HA1E_24H_X1_B23:L19|REP.A025_HA1E_24H_X2_B23:L19|REP.A025_HA1E_24H_X3_B23:L19',
              'pert_id': 'BRD-K77641333',
              'pert_idose': '10.0 um',
              'pert_iname': 'naphazoline',
              'pert_itime': '24 h',
              'pert_type': 'trt_cp',
              'sig_id': 'REP.A025_HA1E_24H:L19'},
 'ncs': {'score': 2.2176711407049905, 'sign': 'negative'},
 'tau': -99.99874925997882,
 'wcs': -0.5047258972246215}

```

To search for the individual signatures belonging to a specific compound (e.g. naphazoline):

``` 
naphazoline_signatures = [item for item in results if item['metadata']['pert_iname'] == 'naphazoline']
``` 

Or, filter the list of results for only signatures profiled in HA1E cells:

```
ha1e_signatures = [item for item in results if item['metadata']['cell_id'] == 'HA1E']
```

Or, filter for signatures with a negative tau score:

```
negative_scores = [item for item in results if item['tau'] < 0.0]
```

Experiments folder: individual comparisons performed in the manuscript