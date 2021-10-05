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

all_disease_UniProtIDs = []
all_disease_IDs = []
all_disease_names = []
all_disease_residues = []
all_disease_changes = []
all_disease_mutation_type = []
disease_changes_combined = []
lenmismatch = []

for path in Paths:
    UniProtID = os.path.basename(os.path.normpath(path))
    
    ## Get disease data
    # =========================================================================
    
    dat_file3 = path+'uprotein.dat'
    
    if os.path.exists(dat_file3):
        text = []
        with open(dat_file3, 'r') as file:
            text = file.read()
            # cut off text when DISEASE_ID starts to begin with NOTE_, Note_ or ClinVarNote_ (as those are only disease notes, not necessarily diseases)
            text = re.split(r'DISEASE_ID.* NOTE_',text,1)[0]
            text = re.split(r'DISEASE_ID.* Note_',text,1)[0]
            text = re.split(r'DISEASE_ID.* ClinVarNote_',text,1)[0]
        
        uniprot_ids = re.findall(r'KEY: (\w+)', text)
        DISEASE_ID = re.findall(r'DISEASE_ID.* (\w+)', text)
        if len(DISEASE_ID)>0:
            DISEASE_NAME = re.findall(r'DISEASE_NAME\[.+\]\[.+\] (\w.*)', text)
            AA_CODE = re.findall(r'AA_CODE.* (.+)', text)
            SEQ_NO = re.findall(r'SEQ_NO.* (\d+)', text)
            AA_MUT = re.findall(r'AA_MUT.* (.+)', text)
            MUT_TYPE = re.findall(r'MUT_TYPE.* (\d)', text)
            
    #        disease = pd.DataFrame(np.column_stack([DISEASE_ID, DISEASE_NAME]), 
    #                                   columns=['DISEASE_ID','DISEASE_NAME'])
    #        disease.to_csv('{}_diseases.tsv'.format(UniProtID), sep='\t', index=False)
            
            if len(AA_MUT)!=len(MUT_TYPE):
                disease_var = pd.DataFrame(np.column_stack([AA_CODE, SEQ_NO, AA_MUT]), 
                                           columns=['AA_CODE','SEQ_NO','AA_MUT'])
                all_disease_mutation_type.append([])
            else:
                disease_var = pd.DataFrame(np.column_stack([AA_CODE, SEQ_NO, AA_MUT, MUT_TYPE]), 
                                           columns=['AA_CODE','SEQ_NO','AA_MUT','MUT_TYPE'])
                all_disease_mutation_type.append(MUT_TYPE)
            
    #        disease_var.to_csv('{}_disease_variants.tsv'.format(UniProtID), sep='\t', index=False)
            if len(AA_MUT)>0:
                disease_changes = disease_var.apply(lambda x: x['AA_CODE']+' -> '+x['AA_MUT'], axis=1)
                disease_changes_list = list(disease_changes)
            else:
                disease_changes_list = []
    #        with open('{}_disease_changes.txt'.format(UniProtID), 'w') as f:
    #            f.write("\n".join(disease_changes_list))
            
            # add to global lists
            all_disease_UniProtIDs.append(UniProtID)
            all_disease_IDs.append(DISEASE_ID)
            all_disease_names.append(DISEASE_NAME)
            all_disease_residues.append(SEQ_NO)
            all_disease_changes.append(disease_changes_list)
            disease_changes_combined.extend(disease_changes_list)

            
## Combine all data and save into files
# =============================================================================

disease_variants_aggregated  = pd.DataFrame(np.column_stack([all_disease_UniProtIDs, all_disease_IDs, all_disease_names, 
                                                             all_disease_residues, all_disease_changes, all_disease_mutation_type]), 
                           columns=['UniProtID','disease_IDs','disease_names','disease_residues','disease_changes','disease_mutation_type'])    

disease_variants_aggregated.to_csv('disease_variants_aggregated.tsv', sep='\t', index=False)

with open('disease_changes_combined.txt', 'w') as f:
    f.write("\n".join(disease_changes_combined))
    
# to open file:
#with open(myfile, 'r') as f:
#    mystring = f.read()
#my_list = mystring.split("\n")

