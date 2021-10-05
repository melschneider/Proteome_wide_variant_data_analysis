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

all_variant_UniProtIDs = []
all_variant_residues = []
all_variant_changes = []
all_variant_max_freq = []
variant_changes_combined = []
lenmismatch = []
novariants = []
all_freqs_df = []

open('treated_uniprots.txt', mode='a').close()

for path in Paths:
    UniProtID = os.path.basename(os.path.normpath(path))
    
    with open("treated_uniprots.txt", "a") as myfile:
        myfile.write(UniProtID+"\n")
    
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
            freq_AFR = re.findall(r'NV_RES_FREQ_AFR\[\d\]\s+(.+)', text)
            freq_AMR = re.findall(r'NV_RES_FREQ_AMR\[\d\]\s+(.+)', text)
            freq_ASJ = re.findall(r'NV_RES_FREQ_ASJ\[\d\]\s+(.+)', text)
            freq_EAS = re.findall(r'NV_RES_FREQ_EAS\[\d\]\s+(.+)', text)
            freq_FIN = re.findall(r'NV_RES_FREQ_FIN\[\d\]\s+(.+)', text)
            freq_NFE = re.findall(r'NV_RES_FREQ_NFE\[\d\]\s+(.+)', text)
            freq_OTH = re.findall(r'NV_RES_FREQ_OTH\[\d\]\s+(.+)', text)
            freq_SAS = re.findall(r'NV_RES_FREQ_SAS\[\d\]\s+(.+)', text)
            
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
            
            if len(variants)==len(freq_AFR)==len(freq_AMR)==len(freq_ASJ)==len(freq_EAS):
                freqs = pd.DataFrame(np.column_stack([variants,synonymous,deleterious,sift_score,polyphen_score,freq_AFR,freq_AMR,freq_ASJ,freq_EAS,freq_FIN,freq_NFE,freq_OTH,freq_SAS]), 
                                               columns=['variants','synonymous','deleterious','SIFT_score','PolyPhen_score','freq_AFR','freq_AMR','freq_ASJ','freq_EAS','freq_FIN','freq_NFE','freq_OTH','freq_SAS'])                
                # remove synonymous variants
                freqs = freqs[freqs['synonymous'] == 'FALSE']
                if len(freqs)>0:
                    # add column with residue number
                    res = freqs['variants'].apply(lambda x: x[3:-3])
                    freqs.insert(loc=1, column='residue', value=res)
                    # add column with residue change
                    changes = freqs['variants'].apply(lambda x: x[0:3]+' -> '+x[-3:])
                    freqs.insert(loc=2, column='change', value=changes)
                    # add column containing UniProtID
                    freqs.insert(loc=0, column='UniProtID', value=UniProtID)
                    # append dataframe to global dataframe
                    all_freqs_df.append(freqs)
                    # save dataframe to file
                    #freqs.to_csv('{}_freqs.tsv'.format(UniProtID), sep='\t', index=False)
                    
                    freqs_only = freqs.loc[:, freqs.columns.str.startswith('freq_')]
                    freqs_only = freqs_only.dropna(how='all')
                    if (len(freqs_only)>0 and isinstance(freqs_only.max().max(), (int, float, complex))):
                        # get the highest frequency of non-synonymous variants, the respective population, mutation, SIFT_score and PolyPhen_score
                        max_freq = [freqs_only.max().max(), 
                                    'Population:'+freqs_only.max().idxmax()[5:], 
                                    'Mutation:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['variants'], 
                                    'SIFT_score:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['SIFT_score'], 
                                    'PolyPhen_score:'+freqs.loc[freqs_only.max(axis=1).idxmax()]['PolyPhen_score']]
                    else:
                        max_freq = []
        
                    # add to global lists
                    all_variant_UniProtIDs.append(UniProtID)
                    all_variant_residues.append(res.to_list())
                    all_variant_changes.append(changes.to_list())
                    all_variant_max_freq.append(max_freq)
                    variant_changes_combined.extend(changes.to_list())
                    
            elif len(variants)==len(synonymous):
                variant_df = pd.DataFrame(np.column_stack([variants,synonymous,deleterious,sift_score,polyphen_score]), 
                                               columns=['variants','synonymous','deleterious','SIFT_score','PolyPhen_score'])                
                # remove synonymous variants
                variant_df = variant_df[variant_df['synonymous'] == 'FALSE']
                if len(variant_df)>0:
                    # add column with residue number
                    res = variant_df['variants'].apply(lambda x: x[3:-3])
                    variant_df.insert(loc=1, column='residue', value=res)
                    # add column with residue change
                    changes = variant_df['variants'].apply(lambda x: x[0:3]+' -> '+x[-3:])
                    variant_df.insert(loc=2, column='change', value=changes)
                    # add column containing UniProtID
                    #variant_df.insert(loc=0, column='UniProtID', value=UniProtID)
                    # save dataframe to file
                    #variant_df.to_csv('{}_variant_df.tsv'.format(UniProtID), sep='\t', index=False)
                    max_freq = []
                    
                    # add to global lists
                    all_variant_UniProtIDs.append(UniProtID)
                    all_variant_residues.append(res.to_list())
                    all_variant_changes.append(changes.to_list())
                    all_variant_max_freq.append(max_freq)
                    variant_changes_combined.extend(changes.to_list())
                    
                    lenmismatch.append(UniProtID+", "+str(len(variants))+", "+str(len(freq_AFR)))
            else:
                novariants.append(UniProtID+", "+str(len(variants))+", "+str(len(variant_df)))
        else:
            novariants.append(UniProtID+", "+str(len(variants)))
        


with open('missing_frequencies.txt', 'w') as f:
    f.write("\n".join(lenmismatch))
with open('no_variants.txt', 'w') as f:
    f.write("\n".join(novariants))

### Combine all data and save into files
## =============================================================================

with open('variant_changes_combined.txt', 'w') as f:
    f.write("\n".join(variant_changes_combined))
    
#variants_aggregated = pd.DataFrame(np.column_stack([all_variant_UniProtIDs, all_variant_residues, all_variant_changes, all_variant_max_freq]), 
#                                   columns=['UniProtID','variant_residues','variant_changes','variant_max_freq'])    

variants_aggregated = pd.DataFrame(np.column_stack([all_variant_UniProtIDs, all_variant_residues, all_variant_changes]), 
                                   columns=['UniProtID','variant_residues','variant_changes'])    

variants_aggregated.to_csv('variants_aggregated.tsv', sep='\t', index=False)

#with open('variant_max_freq.txt', 'w') as f:
#    f.write("\n".join(all_variant_max_freq))

all_freqs_df = pd.concat(all_freqs_df, ignore_index=True)
all_freqs_df.to_csv('variant_freqs_df.tsv', sep='\t', index=False)

# to open file:
#with open(myfile, 'r') as f:
#    mystring = f.read()
#my_list = mystring.split("\n")

#            freq_AFR = re.findall(r'NV_RES_FREQ_AFR\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_AMR = re.findall(r'NV_RES_FREQ_AMR\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_ASJ = re.findall(r'NV_RES_FREQ_ASJ\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_EAS = re.findall(r'NV_RES_FREQ_EAS\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_FIN = re.findall(r'NV_RES_FREQ_FIN\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_NFE = re.findall(r'NV_RES_FREQ_NFE\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_OTH = re.findall(r'NV_RES_FREQ_OTH\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
#            freq_SAS = re.findall(r'NV_RES_FREQ_SAS\[\d\]\s+(\w*\.?\w+[-+]?\w*)', text)
