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
all_variant_UniProtIDs = []
all_variant_residues = []
all_variant_changes = []
all_variant_max_freq = []
variant_changes_combined = []
all_disease_UniProtIDs = []
all_disease_IDs = []
all_disease_names = []
all_disease_residues = []
all_disease_changes = []
all_disease_mutation_type = []
disease_changes_combined = []

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
    
    
    ## Get variant data: vaiant impact (deleterious, SIFT_score, PolyPhen_score) and allele frequencies (gnomAD)
    # =========================================================================
    
    dat_file2 = path+'natvar.dat'
    
    if os.path.exists(dat_file2):
        with open(dat_file2, 'r') as file:
            text = file.read()
        
        variants = re.findall(r'NV_VARIANT\[\d\]\s+(\w+)', text)
        if len(variants)>0:
            synonymous = re.findall(r'NV_SYNONYMOUS\[\d\]\s+(\w+)', text)
            deleterious = re.findall(r'NV_DELETERIOUS\[\d\]\s+(\w+)', text)
            sift_score = re.findall(r'NV_SIFT_SCORE\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            polyphen_score = re.findall(r'NV_POLYPHEN_SCORE\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            #freq_all = re.findall(r'NV_RES_FREQ_ALL\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_AFR = re.findall(r'NV_RES_FREQ_AFR\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_AMR = re.findall(r'NV_RES_FREQ_AMR\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_ASJ = re.findall(r'NV_RES_FREQ_ASJ\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_EAS = re.findall(r'NV_RES_FREQ_EAS\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_FIN = re.findall(r'NV_RES_FREQ_FIN\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_NFE = re.findall(r'NV_RES_FREQ_NFE\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_OTH = re.findall(r'NV_RES_FREQ_OTH\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            freq_SAS = re.findall(r'NV_RES_FREQ_SAS\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
            
            sift_score = [float(x) for x in sift_score if x]
            polyphen_score = [float(x) for x in polyphen_score if x]
            #freq_all = [float(x) for x in freq_all if x]
            freq_AFR = [float(x) for x in freq_AFR if x]
            freq_AMR = [float(x) for x in freq_AMR if x]
            freq_ASJ = [float(x) for x in freq_ASJ if x]
            freq_EAS = [float(x) for x in freq_EAS if x]
            freq_FIN = [float(x) for x in freq_FIN if x]
            freq_NFE = [float(x) for x in freq_NFE if x]
            freq_OTH = [float(x) for x in freq_OTH if x]
            freq_SAS = [float(x) for x in freq_SAS if x]
        
            freqs = pd.DataFrame(np.column_stack([variants,synonymous,deleterious,sift_score,polyphen_score,freq_AFR,freq_AMR,freq_ASJ,freq_EAS,freq_FIN,freq_NFE,freq_OTH,freq_SAS]), 
                                           columns=['variants','synonymous','deleterious','SIFT_score','PolyPhen_score','freq_AFR','freq_AMR','freq_ASJ','freq_EAS','freq_FIN','freq_NFE','freq_OTH','freq_SAS'])
            
            # remove synonymous variants
            freqs = freqs[freqs['synonymous'] == 'FALSE']
            # add column with residue number
            res = freqs['variants'].apply(lambda x: x[3:-3])
            freqs.insert(loc=1, column='residue', value=res)
            # add column with residue change
            changes = freqs['variants'].apply(lambda x: x[0:3]+' -> '+x[-3:])
            freqs.insert(loc=2, column='change', value=changes)
            # add column containing UniProtID
            #freqs.insert(loc=0, column='UniProtID', value=UniProtID)
            # save dataframe to file
            #freqs.to_csv('{}_freqs.tsv'.format(UniProtID), sep='\t', index=False)
            
            freqs_only = freqs.loc[:, freqs.columns.str.startswith('freq_')]
            # get the highest frequency of non-synonymous variants, the respective population, mutation, SIFT_score and PolyPhen_score
            max_freq = [freqs_only.max().max(), 
                        'Population:'+freqs_only.max().idxmax()[5:], 
                        'Mutation:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['variants'], 
                        'SIFT_score:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['SIFT_score'], 
                        'PolyPhen_score:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['PolyPhen_score']]
    
    
            
            # add to global lists
            all_variant_UniProtIDs.append(UniProtID)
            all_variant_residues.append(res.to_list())
            all_variant_changes.append(changes.to_list())
            all_variant_max_freq.append(max_freq)
            variant_changes_combined.extend(changes.to_list())
    

    ## Get disease data
    # =========================================================================
    
    dat_file3 = path+'uprotein.dat'
    
    if os.path.exists(dat_file3):
        text = []
        with open(dat_file3, 'r') as file:
            text = file.read()
            # cut off text when DISEASE_ID starts to begin with NOTE_ or ClinVarNote_ (as those are only disease notes, not necessarily diseases)
            text = re.split(r'DISEASE_ID.* NOTE_',text,1)[0]
            text = re.split(r'DISEASE_ID.* ClinVarNote_',text,1)[0]
        
        uniprot_ids = re.findall(r'KEY: (\w+)', text)
        DISEASE_ID = re.findall(r'DISEASE_ID.* (\w+)', text)
        if len(DISEASE_ID)>0:
            DISEASE_NAME = re.findall(r'DISEASE_NAME\[.+\]\[.+\] (\w.*)', text)
            AA_CODE = re.findall(r'AA_CODE.* (\w+)', text)
            SEQ_NO = re.findall(r'SEQ_NO.* (\d+)', text)
            AA_MUT = re.findall(r'AA_MUT.* (\w+)', text)
            MUT_TYPE = re.findall(r'MUT_TYPE.* (\d)', text)
            
    #        disease = pd.DataFrame(np.column_stack([DISEASE_ID, DISEASE_NAME]), 
    #                                   columns=['DISEASE_ID','DISEASE_NAME'])
    #        disease.to_csv('{}_diseases.tsv'.format(UniProtID), sep='\t', index=False)
    #        
            disease_var = pd.DataFrame(np.column_stack([AA_CODE, SEQ_NO, AA_MUT, MUT_TYPE]), 
                                       columns=['AA_CODE','SEQ_NO','AA_MUT','MUT_TYPE'])
    #        disease_var.to_csv('{}_disease_variants.tsv'.format(UniProtID), sep='\t', index=False)
            
            disease_changes = disease_var.apply(lambda x: x['AA_CODE']+' -> '+x['AA_MUT'], axis=1)
            disease_changes_list = disease_changes.to_list()
            
    #        with open('{}_disease_changes.txt'.format(UniProtID), 'w') as f:
    #            f.write("\n".join(disease_changes_list))
            
            # add to global lists
            all_disease_UniProtIDs.append(UniProtID)
            all_disease_IDs.append(DISEASE_ID)
            all_disease_names.append(DISEASE_NAME)
            all_disease_residues.append(SEQ_NO)
            all_disease_changes.append(disease_changes_list)
            all_disease_mutation_type.append(MUT_TYPE)
            disease_changes_combined.extend(disease_changes_list)


## Combine all data and save into files
# =============================================================================

pdbs_aggregated  = pd.DataFrame(np.column_stack([all_UniProtIDs, all_pdbs]), 
                           columns=['UniProtID','pdbs'])
variants_aggregated  = pd.DataFrame(np.column_stack([all_variant_UniProtIDs, all_variant_residues, all_variant_changes, all_variant_max_freq]), 
                           columns=['UniProtID','variant_residues','variant_changes','variant_max_freq'])    
disease_variants_aggregated  = pd.DataFrame(np.column_stack([all_disease_UniProtIDs, all_disease_IDs, all_disease_names, 
                                                             all_disease_residues, all_disease_changes, all_disease_mutation_type]), 
                           columns=['UniProtID','disease_IDs','disease_names','disease_residues','disease_changes','disease_mutation_type'])    

pdbs_aggregated.to_csv('pdbs_aggregated.tsv', sep='\t', index=False)
variants_aggregated.to_csv('variants_aggregated.tsv', sep='\t', index=False)
disease_variants_aggregated.to_csv('disease_variants_aggregated.tsv', sep='\t', index=False)
with open('variant_changes_combined.txt', 'w') as f:
    f.write("\n".join(variant_changes_combined))
with open('disease_changes_combined.txt', 'w') as f:
    f.write("\n".join(disease_changes_combined))
    


# to open file:
#with open(myfile, 'r') as f:
#    mystring = f.read()
#my_list = mystring.split("\n")

