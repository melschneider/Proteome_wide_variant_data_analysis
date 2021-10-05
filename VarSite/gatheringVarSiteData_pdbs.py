#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:40:07 2021

@author: melanie
"""

import pandas as pd
import numpy as np
import re
import glob
import os.path

Paths = glob.glob('/nfs/research1/thornton/data/DisaStr/uniprot/*/*/')

all_UniProtIDs = []
all_pdbs = []

for path in Paths:
    UniProtID = os.path.basename(os.path.normpath(path))
    
    ## Get PDB-ID and sequence identity for homologous proteins
    # =========================================================================
    
    dat_file1 = path+'pdb.dat'
    if os.path.exists(dat_file1):
        with open(dat_file1, 'r') as file:
            text = file.read()
        
        pdb_ids = re.findall(r'KEY: (\w+)', text)
    #    resolution = re.findall(r'RESOLUTION\s+(\d*\.?\d+)', text)
    #    seq_id_percen = re.findall(r'SEQ_ID_PERCEN\s+(\d*\.?\d+)', text)
    #    seq_id_percen = [float(x) for x in seq_id_percen if x]
    #    e_value = re.findall(r'E_VALUE\s+(\d*\.?\d+e?-?\d+)', text)
    #    e_value = [float(x) for x in e_value if x]
    #
    #    pdbs = pd.DataFrame(np.column_stack([pdb_ids, resolution, seq_id_percen, e_value]), 
    #                                   columns=['pdb_ids','resolution','seq_id_percen','e_value'])
    #    pdbs.to_csv('{}_pdbs.tsv'.format(UniProtID), sep='\t', index=False)
        
        # add to global lists
        all_UniProtIDs.append(UniProtID)
        all_pdbs.append(pdb_ids)
    
    

## Combine all data and save into files
# =============================================================================

pdbs_aggregated  = pd.DataFrame(np.column_stack([all_UniProtIDs, all_pdbs]), 
                           columns=['UniProtID','pdbs'])

pdbs_aggregated.to_csv('pdbs_aggregated.tsv', sep='\t', index=False)

