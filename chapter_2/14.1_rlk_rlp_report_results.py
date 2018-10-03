#!/usr/bin/env python
import pandas as pd
import glob
from functools import reduce
""" 
Python "2.7.6 | 64 -bit"
script is fully functional
looping merge multiple data frames python
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.1'

#RLK
data_rlk = [pd.DataFrame.from_csv(f) for f in glob.glob('*_rlk_modified_results.csv')]
df_final_rlk = reduce(lambda left, right: pd.merge(left, right, how='left', on='domains'), data_rlk)
df_final_transpose_rlk = df_final_rlk.T
df_final_transpose_rlk.to_csv('3_rlk_domain_dataframe.csv')

#RLP
data_rlp = [pd.DataFrame.from_csv(f) for f in glob.glob('*_rlp_modified_results.csv')]
df_final_rlp = reduce(lambda left, right: pd.merge(left, right, how='left', on='domains'), data_rlp)
df_final_transpose_rlp = df_final_rlp.T
df_final_transpose_rlp.to_csv('3_rlp_domain_dataframe.csv')
