#!/usr/bin/env python
import pandas as pd
import glob
""" 
Python "2.7.6 | 64 -bit"
script is fully functional
modified files '*_rlk_total_domains.csv' for join process as a parsing strategy
creates a double column
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.1'

# RLK
column_names_rlk = ['domains', 'domain']
list_of_files_rlk = glob.glob('*_rlk_total_domains.csv')
for filename in list_of_files_rlk:
    df_rlk = pd.read_csv(filename, names=column_names_rlk, header=None)
    df_rlk = df_rlk.assign(domain = df_rlk['domains'])
    df_rlk = df_rlk.rename(columns={'domains': 'domains', 'domain': filename })
    df_rlk = df_rlk.to_csv(filename + '_rlk_modified.csv')
    
# RLP
column_names_rlp = ['domains', 'domain']
list_of_files_rlp = glob.glob('*_rlp_total_domains.csv')
for filename in list_of_files_rlp:
    df_rlp = pd.read_csv(filename, names=column_names_rlk, header=None)
    df_rlp = df_rlp.assign(domain = df_rlp['domains'])
    df_rlp = df_rlp.rename(columns={'domains': 'domains', 'domain': filename })
    df_rlp = df_rlp.to_csv(filename + '_rlp_modified.csv')
