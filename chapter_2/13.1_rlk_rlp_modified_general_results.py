#!/usr/bin/env python
import pandas as pd
""" 
Python "2.7.6 | 64 -bit"
script is fully functional
looping - merge multiple dataframes python
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.1'

df_rlk = pd.read_csv('2_rlk_domain_dataframe.csv')
df_rlk = df_rlk.drop(df_rlk.index[[0]]) # eliminate row 0 on df
final_rlk_df = df_rlk.set_index('domains')
final_rlk_df = final_rlk_df.reindex(["Ins_P5_2-kin", "RIO1", "Pkinase", "PI3_PI4_kinase", 
                                     "Pkinase_Tyr", "Choline_kinase", "ABC1", "Pkinase_C", 
                                     "PIP5K", "WaaY", "APH", "LRRNT_2", "LRR_8", "LRR_1", 
                                     "LRR_4", "LRR_6", "LRR_2", "LRR_5", "Lectin_legB", 
                                     "Lectin_C", "B_lectin", "S_locus_glycop", "PAN_2", 
                                     "PAN_1", "LysM", "Thaumatin", "WAK_assoc", "WAK", 
                                     "GUB_WAK_bind", "Malectin_like", "Malectin", "EGF_CA", 
                                     "EGF", "EGF_3", "Stress-antifung", "RCC1_2", "RVT_2", 
                                     "DUF3403", "Ribonuc_2-5A", "DUF3453", "NAF", "PP2C", 
                                     "PRIMA1", "DUF3660", "GDPD", "rve", "DHHC", "DUF4219", 
                                     "PB1", "Transposase_21", "Transposase_24", "Sec16", 
                                     "gag_pre-integrs", "PI3Ka", "Transpos_assoc", 
                                     "Glyco_hydro_18", "SCAMP", "S1FA", "2OG-FeII_Oxy_2", 
                                     "RETICULATA-like", "MEKHLA", "HPP", "Methyltransf_29", 
                                     "MlaE", "Clathrin", "EDR1", "zf-HIT", "DUF569", "Herpes_gE", 
                                     "Pep3_Vps18", "Beta-lactamase", "Rio2_N", "F-box", "SHQ1", 
                                     "FBD", "MRP-L27", "HSP20", "Terpene_synth_C", "Nodulin_late", 
                                     "Peptidase_M20", "Ost5", "IQ", "Adeno_E3_CR2", "Amino_oxidase", 
                                     "p450", "EF-hand_7", "EMP70", "TMEM154", "OB_NTP_bind", "MatE", 
                                     "EB", "PQQ", "FAD_binding_4", "WD40", "RVP_2", "E1-E2_ATPase", 
                                     "Hydrolase", "EamA", "LEA_2", "DUF4441", "PMR5N", "H_PPase", 
                                     "Cu_bind_like", "Sugar_tr", "Cupin_1", "Retrotran_gag_3", 
                                     "Glyco_hydro_28", "Aldose_epim", "DHQS", "PRA1", "zf-RING_2", 
                                     "Pyrophosphatase", "JmjN", "JmjC", "Peptidase_M50B", "PGG", 
                                     "Ank_2", "SRF-TF", "FATC", "zf-C5HC2", "K-box", "ATP-synt_A", 
                                     "Retrotran_gag_2", "Tmemb_14", "SLAC1"])

final_rlk_df['ATH'] = final_rlk_df['ath_nsg_pfam_parser.csv_rlk_total_domains.csv'].combine_first(final_rlk_df['ath_psg_pfam_parser.csv_rlk_total_domains.csv'])
final_rlk_df.to_csv('A_rlk_domain_dataframe.csv')

# RLP
df_rlp = pd.read_csv('2_rlp_domain_dataframe.csv')
#df_rlp = df_rlp.drop(df_rlp.index[[0]]) # eliminate row 0 on df
final_rlp_df = df_rlp.set_index('domains')
final_rlp_df = final_rlp_df.reindex(["LRR_8", "DUF2854", "LRR_1", "LRR_2", "LRR_4", "LRR_6", "LRR_9", 
                                     "LRR_5", "Gal-bind_lectin", "Glyco_hydro_32C", "XET_C", "Lectin_legB", 
                                     "Glyco_hydro_16", "Calreticulin", "SPRY", "Alginate_lyase2", 
                                     "B_lectin", "S_locus_glycop", "PAN_2", "PAN_1", "PAN_4", "LysM", 
                                     "Thaumatin", "WAK_assoc", "WAK", "GUB_WAK_bind", "Malectin_like", 
                                     "Malectin", "EGF_alliinase", "cEGF", "EGF_CA", "EGF_2", "Stress-antifung", 
                                     "UMP1", "LRRNT_2", "Glyco_hydro_32N", "DUF3357", "DUF604", "Alliinase_C", 
                                     "F-box", "Galactosyl_T", "FBD", "zf-RING_2", "PA", "PRIMA1", "Mito_carr", 
                                     "MtN3_slv", "SHMT", "Peptidase_M8", "Exostosin", "Glyco_transf_90", 
                                     "Ribosomal_S13", "RNA_pol_Rpb2_5", "RNA_pol_Rpb2_4","RNA_pol_Rpb2_7", 
                                     "RNA_pol_Rpb2_6", "RNA_pol_Rpb2_1", "RNA_pol_Rpb2_3","RNA_pol_Rpb2_2", 
                                     "LIAS_N", "DUF4216", "DUF4218", "F-box-like", "Retrotran_gag_2", "Snf7", 
                                     "MatE", "Peptidase_S10", "zf-CCHC", "Glyoxal_oxid_N", "PEARLI-4", "DUF1929"])

final_rlp_df['ATH'] = final_rlp_df['ath_nsg_pfam_parser.csv_rlp_total_domains.csv'].combine_first(final_rlp_df['ath_psg_pfam_parser.csv_rlp_total_domains.csv'])
final_rlp_df.to_csv('A_rlp_domain_dataframe.csv')
