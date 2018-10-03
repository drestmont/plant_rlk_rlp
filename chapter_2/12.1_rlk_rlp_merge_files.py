#!/usr/bin/env python
import pandas as pd
import glob
from functools import reduce
""" 
Python "2.7.6 | 64 -bit"
script is fully functional
merge all '*_rl*_modified.csv' files creating a dataframe using the domain index
as a reference for all species evaluated
looping merge multiple data frames python
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.1'

# RLK
data_rlk = [pd.DataFrame.from_csv(f) for f in glob.glob('*_rlk_modified.csv')]
df_final_rlk = reduce(lambda left, right: pd.merge(left, right, how='left', on='domains'), data_rlk)
df_final_rlk.to_csv('2_rlk_domain_dataframe.csv')

# RLP
data_rlp = [pd.DataFrame.from_csv(f) for f in glob.glob('*_rlp_modified.csv')]
df_final_rlp = reduce(lambda left, right: pd.merge(left, right, how='left', on='domains'), data_rlp)
df_final_rlp.to_csv('2_rlp_domain_dataframe.csv') 
