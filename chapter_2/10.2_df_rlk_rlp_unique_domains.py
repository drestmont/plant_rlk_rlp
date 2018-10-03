#!/usr/bin/env python
import pandas as pd
import glob
""" 
Python "2.7.6 | 64 -bit"
This script is fully functional
read all domains and create a list of unique domains as an index reference
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.2'

results_rlk = pd.DataFrame([])
for counter, file in enumerate(glob.glob('*_rlk_total_domains.csv')):
    name_df_rlk = pd.read_csv(file, header=None)
    results_rlk = results_rlk.append(name_df_rlk)
    
results_rlk.columns = ['domains']
domains_unique_rlk = results_rlk.domains.unique()
domains_unique_list_rlk = domains_unique_rlk.tolist()#total unique groupby domains
domains_unique_join_rlk = ", ".join(domains_unique_list_rlk)# retire external brackets

main_df_rlk = pd.DataFrame({'domains' : domains_unique_list_rlk})
main_df_rlk = main_df_rlk.assign(domain = main_df_rlk['domains'])
main_df_rlk.to_csv('1_rlk_total_unique_domains_rlk_modified.csv')

#RLP
results_rlp = pd.DataFrame([])
for counter, file in enumerate(glob.glob('*_rlp_total_domains.csv')):
    name_df_rlp = pd.read_csv(file, header=None)
    results_rlp = results_rlp.append(name_df_rlp)
    
results_rlp.columns = ['domains']
domains_unique_rlp = results_rlp.domains.unique()
domains_unique_list_rlp = domains_unique_rlp.tolist()#total unique groupby domains
domains_unique_join_rlp = ", ".join(domains_unique_list_rlp)# retire external brackets

main_df_rlp = pd.DataFrame({'domains' : domains_unique_list_rlp})
main_df_rlp = main_df_rlp.assign(domain = main_df_rlp['domains'])
main_df_rlp.to_csv('1_rlp_total_unique_domains_rlp_modified.csv')
