#!/usr/bin/env python
import pandas as pd
import csv
import glob
import re
""" 
Python "2.7.6 | 64 -bit"
This script is fully functional
It parses the Pfam_scan output and identified ectodomains and classified plasma 
membrane proteins as RLKs or RLPs. Each class is independent and combine different
logical operators.

Issues and road map
A. For each plasma membrane protein with its ectodomains it produces a file for
depuration. In the future, it will be producing 1 whole file per species
B. The script follows a sequential approach. Each ectodomain class was built
independently using different combinations. In the future, it will be defined by
functions.
C. Because each outfile makes a file, in the future it will feed a data structure
to build an SQLite db. 

-Input file name
The data include as Input was preprocessed, please read 
methodology: art_chp1_drm_v4_1_editor_reviewed.pdf 
all proteins pre-process have at least 1 transmembrane helix

ath = Arabidopsis thaliana
ath_nsg_pfam_parser.csv
ath_psg_pfam_parser.csv
-ath = species refere
-nsg = no signal peptide / psg=signal peptide.

Output
*_RLK_general_results.csv and *_RLP_general_results.csv produce a summary
The other files report the Protein IDs for each class independently and the sort
of ectodomains present.

ath_nsg_pfam_parser.csv_rlk_general_results.csv = summary
ath_psg_pfam_parser.csv_rlk_general_results.csv = summary

Notes
list of domains reported on pfam reviewed on June 20, 2017
adjusted on January 3, 2018, inclusion of Pkinase_C
    
WARNING
other_d_rlk list must be calculated independently because it depends on the
list of domains identified in the data set evaluated
    
Pending
1. Add example dataset for parsing evaluation
2. Include domains in next version:
LRR: LRRNT(PF01462), LRV(PF01816) & LRRCT(PF01463). Magalhaes et al BMC genomics, 2016
"""
__author__ = "Daniel Restrepo-Montoya"
__version__ = '0.0.9'

################################################################################
# List of ectodomains and cytoplasmic domains
lrr_d = ["DUF285", "FNIP", "LRR_1", "LRR_2", "LRR_3", "LRR_4", "LRR_5", "LRR_6", 
         "LRR_8", "LRR_9", "Recep_L_domain", "LRRNT_2"]

# pkinase clan and Pkinase_C
pkinase_d = ["ABC1", "AceK", "Act-Frag_cataly", "Alpha_kinase", "APH", 
             "APH_6_hur", "Choline_kinase", "CotH", "DUF1679", "DUF2252", 
             "DUF4135", "EcKinase", "Fam20C", "Fructosamin_kin", "FTA2", 
             "Haspin_kinase", "HipA_C", "Ins_P5_2-kin", "IPK", "IucA_IucC",
             "Kdo", "Kinase-like", "KIND", "PI3_PI4_kinase", "PIP49_C", 
             "PIP5K", "Pkinase", "Pkinase_Tyr", "Pox_ser-thr_kin", "RIO1", 
             "Seadorna_VP7", "UL97", "WaaY", "YrbL-PhoP_reg", "YukC", 
             "Pkinase_C"]

l_lectin_d = ["Alginate_lyase2", "ArabFuran-catal", "Bac_rhamnosid", 
              "Calreticulin", "Cleaved_Adhesin", "DUF1080", "DUF1349", 
              "DUF1583", "DUF1961", "DUF2401", "DUF3472", "DUF4975", 
              "Exotox-A_bind", "Gal-bind_lectin", "Glyco_hydro_11", 
              "Glyco_hydro_12", "Glyco_hydro_16", "Glyco_hydro_32C", 
              "Glyco_hydro_7", "Laminin_G_1", "Laminin_G_2", "Laminin_G_3", 
              "Lectin_leg-like", "Lectin_legB", "MAM", "Methyltransf_FA", 
              "Neuralized", "Pentaxin", "Peptidase_A4", "Polysacc_lyase", 
              "PRY", "Reoviridae_Vp9", "Sial-lect-inser", "Sialidase", "SKN1", 
              "Spike_NTD", "SPRY", "TgMIC1", "Toxin_R_bind_N", "TSP_C", 
              "VP4_haemagglut", "XET_C", "YrpD"]
              
c_lectin_d = ["Lectin_C"]

g_lectin_b_d = ["B_lectin"]

g_lectin_s_d = ["S_locus_glycop"]

g_lectin_p_d = ["AMA-1", "MANEC", "PAN_1", "PAN_2", "PAN_3", "PAN_4"]

lysm_d = ["LysM", "OapA", "Phage_tail_X"]

pr5k_d = ["Thaumatin"]

tnfr_d = ["BaffR-Tall_bind", "BCMA-Tall_bind", "NCD3G", "stn_TNFRSF12A", 
          "TACI-CRD2", "TNFR_c6"]
           
wak_d = ["WAK", "GUB_WAK_bind", "WAK_assoc"]

malectin_d = ["Malectin", "Malectin_like"]

egf_d = ["cEGF", "CFC", "DSL", "EGF", "EGF_2", "EGF_3", "EGF_alliinase", 
         "EGF_CA", "EGF_MSP1_1", "FOLN", "FXa_inhibition", "Gla", "hEGF", 
         "Laminin_EGF", "Plasmod_Pvs28", "Sushi", "Sushi_2", "Tme5_EGF_like"]

stress_antifung_d = ["Stress-antifung"]

all_d = ["DUF285", "FNIP", "LRR_1", "LRR_2", "LRR_3", "LRR_4", "LRR_5", 
         "LRR_6", "LRR_8", "LRR_9", "Recep_L_domain", "LRRNT_2", "ABC1", 
         "AceK", "Act-Frag_cataly", "Alpha_kinase", "APH", "APH_6_hur", 
         "Choline_kinase", "CotH", "DUF1679", "DUF2252", "DUF4135", "EcKinase",
         "Fam20C", "Fructosamin_kin", "FTA2", "Haspin_kinase", "HipA_C", 
         "Ins_P5_2-kin", "IPK", "IucA_IucC", "Kdo", "Kinase-like", "KIND", 
         "PI3_PI4_kinase", "PIP49_C", "PIP5K", "Pkinase", "Pkinase_Tyr", 
         "Pox_ser-thr_kin", "RIO1", "Seadorna_VP7", "UL97", "WaaY", 
         "YrbL-PhoP_reg", "YukC", "Alginate_lyase2", "ArabFuran-catal", 
         "Bac_rhamnosid", "Calreticulin", "Cleaved_Adhesin", "DUF1080", 
         "DUF1349", "DUF1583", "DUF1961", "DUF2401", "DUF3472", "DUF4975", 
         "Exotox-A_bind", "Gal-bind_lectin", "Glyco_hydro_11", "Glyco_hydro_12", 
         "Glyco_hydro_16", "Glyco_hydro_32C", "Glyco_hydro_7", "Laminin_G_1", 
         "Laminin_G_2", "Laminin_G_3", "Lectin_leg-like", "Lectin_legB", "MAM", 
         "Methyltransf_FA", "Neuralized", "Pentaxin", "Peptidase_A4", 
         "Polysacc_lyase", "PRY", "Reoviridae_Vp9", "Sial-lect-inser", 
         "Sialidase", "SKN1", "Spike_NTD", "SPRY", "TgMIC1", "Toxin_R_bind_N", 
         "TSP_C", "VP4_haemagglut", "XET_C", "YrpD", "Lectin_C", "B_lectin", 
         "S_locus_glycop", "AMA-1", "MANEC", "PAN_1", "PAN_2", "PAN_3", "PAN_4", 
         "LysM", "OapA", "Phage_tail_X", "Thaumatin", "BaffR-Tall_bind", 
         "BCMA-Tall_bind", "NCD3G", "stn_TNFRSF12A", "TACI-CRD2", "TNFR_c6", 
         "WAK", "GUB_WAK_bind", "WAK_assoc","Malectin", "Malectin_like", "cEGF", 
         "CFC", "DSL", "EGF", "EGF_2", "EGF_3", "EGF_alliinase", "EGF_CA", 
         "EGF_MSP1_1", "FOLN", "FXa_inhibition", "Gla", "hEGF", "Laminin_EGF", 
         "Plasmod_Pvs28", "Sushi", "Sushi_2", "Tme5_EGF_like", "Stress-antifung"]

# no pkinase domains included - all ectodomain target (LRR, G/L/C-lectin, LysM, pr5k,
# Malectin, Wak, EGF, TNFR, Stress-antifung
almost_all_d = ["DUF285", "FNIP", "LRR_1", "LRR_2", "LRR_3", "LRR_4", "LRR_5", 
                "LRR_6", "LRR_8", "LRR_9", "Recep_L_domain", "LRRNT_2", 
                "Alginate_lyase2", "ArabFuran-catal", "Bac_rhamnosid", 
                "Calreticulin", "Cleaved_Adhesin", "DUF1080", "DUF1349", 
                "DUF1583", "DUF1961", "DUF2401", "DUF3472", "DUF4975", 
                "Exotox-A_bind", "Gal-bind_lectin", "Glyco_hydro_11", 
                "Glyco_hydro_12", "Glyco_hydro_16", "Glyco_hydro_32C", 
                "Glyco_hydro_7", "Laminin_G_1", "Laminin_G_2", "Laminin_G_3", 
                "Lectin_leg-like", "Lectin_legB", "MAM", "Methyltransf_FA", 
                "Neuralized", "Pentaxin", "Peptidase_A4", "Polysacc_lyase", 
                "PRY", "Reoviridae_Vp9", "Sial-lect-inser", "Sialidase", "SKN1", 
                "Spike_NTD", "SPRY", "TgMIC1", "Toxin_R_bind_N", "TSP_C", 
                "VP4_haemagglut", "XET_C", "YrpD", "Lectin_C", "B_lectin", 
                "S_locus_glycop", "AMA-1", "MANEC", "PAN_1", "PAN_2", "PAN_3", 
                "PAN_4", "LysM", "OapA", "Phage_tail_X", "Thaumatin", 
                "BaffR-Tall_bind", "BCMA-Tall_bind", "NCD3G", "stn_TNFRSF12A", 
                "TACI-CRD2", "TNFR_c6", "WAK", "GUB_WAK_bind", "WAK_assoc", 
                "Malectin", "Malectin_like", "cEGF", "CFC", "DSL", "EGF", 
                "EGF_2", "EGF_3", "EGF_alliinase", "EGF_CA", "EGF_MSP1_1", 
                "FOLN", "FXa_inhibition", "Gla", "hEGF", "Laminin_EGF", 
                "Plasmod_Pvs28", "Sushi", "Sushi_2", "Tme5_EGF_like", 
                "Stress-antifung"]

# NB-ARC http://pfam.xfam.org/clan/CL0023
nb_arc_d = ["6PF2K", "AAA", "AAA-ATPase_like", "AAA_10", "AAA_11", "AAA_12", 
            "AAA_13", "AAA_14", "AAA_15", "AAA_16", "AAA_17", "AAA_18", 
            "AAA_19", "AAA_2", "AAA_21", "AAA_22", "AAA_23", "AAA_24", 
            "AAA_25", "AAA_26", "AAA_27", "AAA_28", "AAA_29", "AAA_3", "AAA_30", 
            "AAA_31", "AAA_32", "AAA_33", "AAA_34", "AAA_35", "AAA_5", "AAA_6", 
            "AAA_7", "AAA_8", "AAA_9", "AAA_PrkA", "ABC_ATPase", "ABC_tran", 
            "ABC_tran_Xtn", "Adeno_IVa2", "Adenylsucc_synt", "ADK", 
            "AFG1_ATPase", "AIG1", "APS_kinase", "Arf", "ArgK", "ArsA_ATPase", 
            "ATP-synt_ab", "ATP_bind_1", "ATP_bind_2", "ATPase", "ATPase_2", 
            "Bac_DnaA", "BCA_ABC_TP_C", "Beta-Casp", "Cas_Csn2", "Cas_St_Csn2", 
            "CbiA", "CBP_BcsQ", "CDC73_C", "CENP-M", "CFTR_R", "CLP1_P", "CMS1", 
            "CoaE", "CobA_CobO_BtuR", "CobU", "cobW", "CPT", "CSM2", 
            "CTP_synth_N", "Cytidylate_kin", "Cytidylate_kin2", "DAP3", "DEAD", 
            "DEAD_2", "DLIC", "DNA_pack_C", "DNA_pack_N", "DNA_pol3_delta", 
            "DNA_pol3_delta2", "DnaB_C", "dNK", "DUF1611", "DUF1726", 
            "DUF2075", "DUF2326", "DUF2478", "DUF257", "DUF2791", "DUF2813", 
            "DUF3584", "DUF463", "DUF815", "DUF853", "DUF87", "DUF927", 
            "Dynamin_N", "Dynein_heavy", "ERCC3_RAD25_C", "Exonuc_V_gamma", 
            "FeoB_N", "Fer4_NifH", "Flavi_DEAD", "FTHFS", "FtsK_SpoIIIE", 
            "G-alpha", "Gal-3-0_sulfotr", "GBP", "GBP_C", "GTP_EFTU", 
            "Gtr1_RagA", "Guanylate_kin", "GvpD", "HDA2-3", "Helicase_C", 
            "Helicase_C_2", "Helicase_C_4", "Helicase_RecD", "Herpes_Helicase", 
            "Herpes_ori_bp", "Herpes_TK", "Hydin_ADK", "IIGP", "IPPT", "IPT", 
            "IstB_IS21", "KAP_NTPase", "KdpD", "Kinesin", "KTI12", "LAP1C", 
            "Lon_2", "LpxK", "MCM", "MEDS", "Mg_chelatase", "Microtub_bd", 
            "MipZ", "MMR_HSR1", "MMR_HSR1_C", "MobB", "MukB", "MutS_V", 
            "Myosin_head", "NACHT", "NB-ARC", "NOG1", "NTPase_1", "NTPase_P4", 
            "ORC3_N", "ParA", "Parvo_NS1", "PAXNEB", "PduV-EutP", "PhoH", 
            "PIF1", "Podovirus_Gp16", "Polyoma_lg_T_C", "Pox_A32", "PPK2", 
            "PPV_E1_C", "PRK", "PSY3", "Rad17", "Rad51", "Ras", "RecA", 
            "ResIII", "RHD3", "RHSP", "RNA12", "RNA_helicase", "Roc", 
            "RsgA_GTPase", "RuvB_N", "SbcCD_C", "SecA_DEAD", "Septin", 
            "Sigma54_activ_2", "Sigma54_activat", "SKI", "SMC_N", "SNF2_N", 
            "Spore_IV_A", "SRP54", "SRPRB", "SulA", "Sulfotransfer_1", 
            "Sulfotransfer_2", "Sulfotransfer_3", "Sulphotransf", "T2SSE", 
            "T4SS-DNA_transf", "Terminase_1", "Terminase_3", "Terminase_6", 
            "Terminase_GpA", "Thymidylate_kin", "TIP49", "TK", "TniB", "Torsin", 
            "TraG-D_C", "tRNA_lig_kinase", "TrwB_AAD_bind", "TsaE", "UvrB", 
            "UvrD-helicase", "UvrD_C", "UvrD_C_2", "Viral_helicase1", "VirC1", 
            "VirE", "Zeta_toxin", "Zot"]

# after identification of all domains present in all datasets rlk and rlp, exclude 
# target ectodomains                        
other_d_rlk = ["Aph-1", "RCC1_2", "DUF3403", "Ribonuc_2-5A", "NAF", "DUF3660", 
               "Glyco_hydro_18", "RVT_2", "DUF3453", "PP2C", "PRIMA1", "GDPD", 
               "rve", "DHHC", "DUF4219", "PB1", "Transposase_21", 
               "Transposase_24", "Sec16", "gag_pre-integrs", "PI3Ka", 
               "Transpos_assoc", "SCAMP", "S1FA", "2OG-FeII_Oxy_2", 
               "RETICULATA-like", "MEKHLA", "HPP", "Methyltransf_29", 
               "MlaE", "Clathrin", "EDR1", "zf-HIT", "DUF569", "Herpes_gE", 
               "Pep3_Vps18", "Beta-lactamase", "Rio2_N", "F-box", "SHQ1", "FBD", 
               "MRP-L27", "HSP20", "Terpene_synth_C", "Nodulin_late", 
               "Peptidase_M20", "Ost5", "IQ", "Adeno_E3_CR2", "Amino_oxidase", 
               "p450", "EF-hand_7", "EMP70", "TMEM154", "OB_NTP_bind", "MatE", 
               "EB", "PQQ", "FAD_binding_4", "WD40", "RVP_2", "E1-E2_ATPase", 
               "Hydrolase", "EamA", "LEA_2", "DUF4441", "PMR5N", "H_PPase", 
               "Cu_bind_like", "Sugar_tr", "Cupin_1", "Retrotran_gag_3", 
               "Glyco_hydro_28", "Aldose_epim", "DHQS", "PRA1", "zf-RING_2", 
               "Pyrophosphatase", "JmjN", "JmjC", "Peptidase_M50B", "PGG", 
               "Ank_2", "SRF-TF", "FATC", "zf-C5HC2", "K-box", "ATP-synt_A", 
               "Retrotran_gag_2", "Tmemb_14", "SLAC1"]

# process RLP and RLK - convert each row in the file to unique id and a vector 
# of domains             
list_of_files = glob.glob('*.csv')
for filename in list_of_files:
    df = pd.read_csv(filename)
    df_group_id = df.groupby(('seq_id'))
    res_rlk_rlp = df_group_id['hmm_name'].unique()
    res_rlk_rlp = res_rlk_rlp.reset_index()
    res_rlk_rlp = res_rlk_rlp.set_index('seq_id')
    res_hmm_name = res_rlk_rlp['hmm_name'].astype(str)
    
    ################# RLP filters ###################

    # RLP_LRR_ NOT pkinase NOT nb-ARC
    # evaluate true/false condition
    set_rlp_lrr = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lrr_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    # subset true condition - convert to string
    set_rlp_lrr = (res_hmm_name[res_hmm_name.str.contains('|'.join(lrr_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    # extract id from index - vector
    set_idlist_rlp_lrr = set_rlp_lrr.index.tolist()

    with open(filename + "_rlp_lrr_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_lrr)
    
    df_rlp_lrr = pd.read_csv(filename + "_rlp_lrr_df.csv")
    df_rlp_lrr["lrr"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_lrr_df.csv", index=False)


    # RLP_l-lectin: l-lectin NOT pkinase NOT nb-ARC
    set_rlp_llectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(l_lectin_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_llectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(l_lectin_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_llectin = set_rlp_llectin.index.tolist()

    with open (filename + "_rlp_l-lectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_llectin)
    
    df_rlp_lrr = pd.read_csv(filename + "_rlp_l-lectin_df.csv")
    df_rlp_lrr["l-lectin"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_l-lectin_df.csv", index=False)
    
    
    # RLP_c-lectin: c-lectin NOT pkinase NOT nb-ARC
    set_rlp_clectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(c_lectin_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_clectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(c_lectin_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_clectin = set_rlp_clectin.index.tolist()

    with open (filename + "_rlp_c-lectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_clectin)
        
    df_rlp_lrr = pd.read_csv(filename + "_rlp_c-lectin_df.csv")
    df_rlp_lrr["c-lectin"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_c-lectin_df.csv", index=False)    
    
    ################################################
    # RLP_g-lectin-b: g-lectin-b NOT pkinase NOT nb-ARC
    set_rlp_glectin_b = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_b_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_glectin_b = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_b_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_glectin_b = set_rlp_glectin_b.index.tolist()

    with open (filename + "_rlp_g-lectin-b_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_glectin_b)
   
    df_rlp_lrr = pd.read_csv(filename + "_rlp_g-lectin-b_df.csv")
    df_rlp_lrr["g-lectin-b"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_g-lectin-b_df.csv", index=False)
        
                  
    # RLP_g-lectin-s: g-lectin-s NOT pkinase NOT nb-ARC
    set_rlp_glectin_s = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_s_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_glectin_s = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_s_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_glectin_s = set_rlp_glectin_s.index.tolist()

    with open (filename + "_rlp_g-lectin-s_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_glectin_s)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_g-lectin-s_df.csv")
    df_rlp_lrr["g-lectin-s"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_g-lectin-s_df.csv", index=False)
        
    # RLP_g-lectin-p: g-lectin-p NOT pkinase NOT nb-ARC
    set_rlp_glectin_p = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_p_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_glectin_p = (res_hmm_name[res_hmm_name.str.contains('|'.join(g_lectin_p_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_glectin_p = set_rlp_glectin_p.index.tolist()

    with open (filename + "_rlp_g-lectin-p_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_glectin_p)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_g-lectin-p_df.csv")
    df_rlp_lrr["g-lectin-p"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_g-lectin-p_df.csv", index=False)
       
    #################################################    
    # RLP_lysm: lysm NOT pkinase NOT nb-ARC
    set_rlp_lysm = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lysm_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_lysm = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lysm_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_lysm = set_rlp_lysm.index.tolist()

    with open (filename + "_rlp_lysm_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_lysm)
    
    df_rlp_lrr = pd.read_csv(filename + "_rlp_lysm_df.csv")
    df_rlp_lrr["lysm"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_lysm_df.csv", index=False)
    
    # RLP_pr5k: pr5k NOT pkinase NOT nb-ARC
    set_rlp_pr5k = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pr5k_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_pr5k = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pr5k_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_pr5k = set_rlp_pr5k.index.tolist()

    with open (filename + "_rlp_pr5k_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_pr5k)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_pr5k_df.csv")
    df_rlp_lrr["pr5k"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_pr5k_df.csv", index=False)
        
    # RLP_tnfr: tnfr NOT pkinase NOT nb-ARC
    set_rlp_tnfr = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(tnfr_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_tnfr = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(tnfr_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_tnfr = set_rlp_tnfr.index.tolist()

    with open (filename + "_rlp_tnfr_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_tnfr)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_tnfr_df.csv")
    df_rlp_lrr["tnfr"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_tnfr_df.csv", index=False)
    
    # RLP_wak: wak NOT pkinase NOT nb-ARC
    set_rlp_wak = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(wak_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_wak = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(wak_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_wak = set_rlp_wak.index.tolist()

    with open (filename + "_rlp_wak_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_wak)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_wak_df.csv")
    df_rlp_lrr["wak"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_wak_df.csv", index=False)
        
    # RLP_malectin: malectin NOT pkinase NOT nb-ARC
    set_rlp_malectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(malectin_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_malectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(malectin_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_malectin = set_rlp_malectin.index.tolist()

    with open (filename + "_rlp_malectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_malectin)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_malectin_df.csv")
    df_rlp_lrr["malectin"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_malectin_df.csv", index=False)
        
    # RLP_egf: egf NOT pkinase NOT nb-ARC
    set_rlp_egf = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(egf_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_egf = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(egf_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_egf = set_rlp_egf.index.tolist()

    with open (filename + "_rlp_egf_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_egf)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_egf_df.csv")
    df_rlp_lrr["egf"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_egf_df.csv", index=False)
        
    # RLP_egf: stress_antifung NOT pkinase NOT nb-ARC
    set_rlp_str_anti = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(stress_antifung_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlp_str_anti = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(stress_antifung_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlp_str_anti = set_rlp_str_anti.index.tolist()

    with open (filename + "_rlp_stress_antifung_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlp_str_anti)

    df_rlp_lrr = pd.read_csv(filename + "_rlp_stress_antifung_df.csv")
    df_rlp_lrr["stress_antifung"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlp_stress_antifung_df.csv", index=False)
        
    # concatenate all filters - identify all domains present on each dataset per species
    frames_rlp = [set_rlp_lrr, set_rlp_llectin, set_rlp_clectin, set_rlp_glectin_b, \
                  set_rlp_glectin_s, set_rlp_glectin_p, set_rlp_lysm, set_rlp_pr5k, \
                  set_rlp_tnfr, set_rlp_wak, set_rlp_malectin, set_rlp_egf, \
                  set_rlp_str_anti]
    
    total_rlp = pd.concat(frames_rlp, verify_integrity=False) # True integrity report duplicates
    total_rlp_unique = total_rlp.unique() # set of domains identified
    total_rlp_unique_list = total_rlp_unique.tolist()#total unique groupby domains
    total_rlp_join = ", ".join(total_rlp_unique_list)# retire external brackets
    total_list_rlp = re.findall(r"'?\w[\w']*(?:-\w+)*'?'", total_rlp_join)
    total_list_rlp_mod = [item.replace("'", "") for item in total_list_rlp]
    total_list_rlp_unique = list(set(total_list_rlp_mod))#items unique
    #total_list_words = '\n'.join(total_list_unique)
    
    with open (filename + "_rlp_total_domains.csv", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(total_list_rlp_unique)
        f.close()

    # create a summary file with results per class - results could be redundant use it as reference     
    with open (filename + "_rlp_general_results.csv", "w") as f:
        rlp_lrr = set_rlp_lrr.count()
        f.write( 'rlp_lrr,' + str(rlp_lrr) + '\n')
        
        rlp_llectin = set_rlp_llectin.count()
        f.write( 'rlp_l-lectine,' + str(rlp_llectin) + '\n')
        
        rlp_clectin = set_rlp_clectin.count()
        f.write( 'rlp_c-lectin,' + str(rlp_clectin) + '\n')
        
        rlp_glectin_b = set_rlp_glectin_b.count()
        f.write( 'rlp_g-lectin-b,' + str(rlp_glectin_b) + '\n')
        
        rlp_glectin_s = set_rlp_glectin_s.count()
        f.write( 'rlp_g-lectin-s,' + str(rlp_glectin_s) + '\n')
        
        rlp_glectin_p = set_rlp_glectin_p.count()
        f.write( 'rlp_g-lectin-p,' + str(rlp_glectin_p) + '\n')
        
        rlp_lysm = set_rlp_lysm.count()
        f.write('rlp_lysm,' + str(rlp_lysm) + '\n')
        
        rlp_malectin = set_rlp_malectin.count()
        f.write('rlp_malectin,' + str(rlp_malectin) + '\n')
        
        rlp_pr5k = set_rlp_pr5k.count()
        f.write('rlp_pr5k,' + str(rlp_pr5k) + '\n')
        
        rlp_tnfr = set_rlp_tnfr.count()
        f.write('rlp_tnfr,' + str(rlp_tnfr) + '\n')
        
        rlp_wak = set_rlp_wak.count()
        f.write('rlp_wak,' + str(rlp_wak) + '\n')
        
        rlp_egf = set_rlp_egf.count()
        f.write('rlp_egf,' + str(rlp_egf) + '\n')
        
        rlp_stress_antifug = set_rlp_str_anti.count()
        f.write('rlp_stress_antifug,' + str(rlp_stress_antifug) + '\n')

        rlp_total = (rlp_lrr + rlp_llectin + rlp_clectin + rlp_glectin_b 
                     + rlp_glectin_s + rlp_glectin_p + rlp_lysm + rlp_pr5k 
                     + rlp_wak + rlp_malectin + rlp_egf + rlp_stress_antifug)
        f.write('rlp_total,' + str(rlp_total) + '\n')

    ################################# RLK ######################################
    # rlk_LRR: lrr AND pkinase
    # evaluate true/false condition
    set_rlk_lrr = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lrr_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    # subset true condition - convert to string
    set_rlk_lrr = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lrr_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    #set_rlk_lrr.to_csv(filename + '_rlk_lrr_df.csv')
    # extract id from index - vector
    set_idlist_rlk_lrr = set_rlk_lrr.index.tolist()

    with open(filename + "_rlk_lrr_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_lrr)
 
    df_rlp_lrr = pd.read_csv(filename + "_rlk_lrr_df.csv")
    df_rlp_lrr["lrr"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlk_lrr_df.csv", index=False)
   
    
    # rlk_l-lectin: l-lectin AND pkinase
    set_rlk_llectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(l_lectin_d)) 
                       & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_llectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(l_lectin_d)) 
                       & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_llectin = set_rlk_llectin.index.tolist()

    with open (filename + "_rlk_l-lectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_llectin)
    
    df_rlp_lrr = pd.read_csv(filename + "_rlk_l-lectin_df.csv")
    df_rlp_lrr["l-lectin"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlk_l-lectin_df.csv", index=False)
    
    # rlk_c-lectin: c-lectin AND pkinase
    set_rlk_clectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(c_lectin_d)) 
                       & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_clectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(c_lectin_d)) 
                       & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                       & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_clectin = set_rlk_clectin.index.tolist()

    with open (filename + "_rlk_c-lectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_clectin)
    
    df_rlp_lrr = pd.read_csv(filename + "_rlk_c-lectin_df.csv")
    df_rlp_lrr["c-lectin"] = "1"
    df_rlp_lrr.to_csv(filename + "_rlk_c-lectin_df.csv", index=False)
    
    ##############################################################       
    # rlk_g-lectin_b: g-lectin b-lectin AND pkinase
    set_rlk_glectin_b = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_b_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_glectin_b = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_b_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_glectin_b = set_rlk_glectin_b.index.tolist()

    with open (filename + "_rlk_g-lectin-b_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_glectin_b)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_g-lectin-b_df.csv")
    df_rlk_lrr["g-lectin-b"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_g-lectin-b_df.csv", index=False)
    
    # rlk_g-lectin_s: g-lectin s-locus AND pkinase
    set_rlk_glectin_s = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_s_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_glectin_s = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_s_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_glectin_s = set_rlk_glectin_s.index.tolist()

    with open (filename + "_rlk_g-lectin-s_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_glectin_s)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_g-lectin-s_df.csv")
    df_rlk_lrr["g-lectin-s"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_g-lectin-s_df.csv", index=False)
   
     # rlk_g-lectin_b:  PAN domain AND pkinase
    set_rlk_glectin_p = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_p_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_glectin_p = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(g_lectin_p_d)) 
                         & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                         & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_glectin_p = set_rlk_glectin_p.index.tolist()

    with open (filename + "_rlk_g-lectin-p_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_glectin_p)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_g-lectin-p_df.csv")
    df_rlk_lrr["g-lectin-p"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_g-lectin-p_df.csv", index=False)
    
    ############################################################################    
    # rlk_lysm: lysm AND pkinase
    set_rlk_lysm = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lysm_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_lysm = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(lysm_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_lysm = set_rlk_lysm.index.tolist()

    with open (filename + "_rlk_lysm_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_lysm)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_lysm_df.csv")
    df_rlk_lrr["lysm"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_lysm_df.csv", index=False)
        
    # rlk_pr5k: pr5k AND pkinase
    set_rlk_pr5k = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pr5k_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_pr5k = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pr5k_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_pr5k = set_rlk_pr5k.index.tolist()

    with open (filename + "_rlk_pr5k_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_pr5k)
   
    df_rlk_lrr = pd.read_csv(filename + "_rlk_pr5k_df.csv")
    df_rlk_lrr["pr5k"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_pr5k_df.csv", index=False)
        
    # rlk_tnfr: tnfr AND pkinase
    set_rlk_tnfr = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(tnfr_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_tnfr = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(tnfr_d)) 
                    & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_tnfr = set_rlk_tnfr.index.tolist()

    with open (filename + "_rlk_tnfr_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_tnfr)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_tnfr_df.csv")
    df_rlk_lrr["tnfr"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_tnfr_df.csv", index=False)
    
    # rlk_wak: wak AND pkinase
    set_rlk_wak = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(wak_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_wak = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(wak_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_wak = set_rlk_wak.index.tolist()

    with open (filename + "_rlk_wak_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_wak)

    df_rlk_lrr = pd.read_csv(filename + "_rlk_wak_df.csv")
    df_rlk_lrr["wak"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_wak_df.csv", index=False)
        
    # rlk_malectin: malectin AND pkinase
    set_rlk_malectin = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(malectin_d)) 
                        & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_malectin = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(malectin_d)) 
                        & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_malectin = set_rlk_malectin.index.tolist()

    with open (filename + "_rlk_malectin_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_malectin)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_malectin_df.csv")
    df_rlk_lrr["malectin"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_malectin_df.csv", index=False)
        
    # rlk_egf: egf AND pkinase
    set_rlk_egf = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(egf_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_egf = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(egf_d)) 
                   & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                   & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_egf = set_rlk_egf.index.tolist()

    with open (filename + "_rlk_egf_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_egf)
    
    df_rlk_lrr = pd.read_csv(filename + "_rlk_egf_df.csv")
    df_rlk_lrr["egf"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_egf_df.csv", index=False)
        
    # rlk_egf: stress_antifung NOT pkinase
    set_rlk_str_anti = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(stress_antifung_d)) 
                        & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)))
    set_rlk_str_anti = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(stress_antifung_d)) 
                        & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                        & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    set_idlist_rlk_str_anti = set_rlk_str_anti.index.tolist()

    with open (filename + "_rlk_stress_antifung_df.csv", "w") as f:
        wr = csv.writer(f, delimiter="\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_str_anti)
        
    df_rlk_lrr = pd.read_csv(filename + "_rlk_stress_antifung_df.csv")
    df_rlk_lrr["stress_antifung"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_stress_antifung_df.csv", index=False)
        
    ############################################################################
    # rlk_extra pkinase positive AND signal peptide NOT ectodomain target NOT nb-arc= other RLK
    # rlk_extra pkinase negative signal peptide NOT ectodomain target NOT nb-arc = probable RLCK
    # not lrr, l_lectin, c_lectin, g_lectin, lysm, pr5k, tnfr, wak, malectin, egf, stress_antifung, nb-arc
    set_rlk_extra_or_rlck = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                             & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(other_d_rlk)) 
                             & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(almost_all_d)) 
                             & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))) 
    #set_rlk_extra_or_rlck.to_csv('1_proof')
    set_rlk_extra_or_rlck = (res_hmm_name[res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                             & res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(other_d_rlk)) 
                             & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(almost_all_d)) 
                             & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d))])
    #set_rlk_extra.to_csv(filename + '_rlk_extra_df.csv')
    set_idlist_rlk_extra_or_rlck = set_rlk_extra_or_rlck.index.tolist()

    with open (filename + "_rlk_extra_or_rlck_df.csv", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_extra_or_rlck)
        
    df_rlk_lrr = pd.read_csv(filename + "_rlk_extra_or_rlck_df.csv")
    df_rlk_lrr["rlk-extra-or-rlck"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_extra_or_rlck_df.csv", index=False)
        
    # identify unique domains present in rlk_kinases_other_domains
    ref_rlk_extra_or_rlck = set_rlk_extra_or_rlck.unique()
    ref_rlk_extra_or_rlck_list = ref_rlk_extra_or_rlck.tolist()
    
    with open (filename + "_refdom_rlk_extra_or_rlck.csv", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(["protein_id"])
        wr.writerow(ref_rlk_extra_or_rlck_list)
        f.close()
    
    ############################################################################
    # RLCK: pkinase NOT signal peptide NOT ectodomain target NOT nb-arc NOT other domains different to pkinase
    set_rlk_or_rlck_only_pkinase = (res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(pkinase_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(almost_all_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(other_d_rlk)))
    set_rlk_or_rlck_only_pkinase = (res_hmm_name[res_hmm_name.str.contains('|'.join(pkinase_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(almost_all_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(nb_arc_d)) 
                                    & ~res_hmm_name.str.contains(r'(?<!-)\b(?:%s)\b(?!-)' % '|'.join(other_d_rlk))])
    #set_rlk_extra.to_csv(filename + '_rlk_extra_df.csv')
    set_idlist_rlk_or_rlck_only_pkinase = set_rlk_or_rlck_only_pkinase.index.tolist()

    with open (filename + "_rlk_or_rlck_only_pkinase_df.csv", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(["protein_id"])
        wr.writerow(set_idlist_rlk_or_rlck_only_pkinase)
        
    df_rlk_lrr = pd.read_csv(filename + "_rlk_or_rlck_only_pkinase_df.csv")
    df_rlk_lrr["rlk_or_rlck_only_pkinase"] = "1"
    df_rlk_lrr.to_csv(filename + "_rlk_or_rlck_only_pkinase_df.csv", index=False)
    
    # concatenate all filters
    frames_rlk = [set_rlk_lrr, set_rlk_llectin, set_rlk_clectin, set_rlk_glectin_b, \
                  set_rlk_glectin_s, set_rlk_glectin_p, set_rlk_lysm, set_rlk_pr5k, \
                  set_rlk_tnfr, set_rlk_wak, set_rlk_malectin, set_rlk_egf, \
                  set_rlk_str_anti, set_rlk_extra_or_rlck, set_rlk_or_rlck_only_pkinase]
    
    total_rlk = pd.concat(frames_rlk, verify_integrity=False) # True integrity report duplicates
    total_rlk_unique = total_rlk.unique() # set of domains identified
    total_rlk_unique_list = total_rlk_unique.tolist()#total unique groupby domains
    total_rlk_join = ", ".join(total_rlk_unique_list)# retire external brackets
    #total_list_rlk = re.findall(r"'(\w+)'", total_rlk_join)#list only one brackets
    total_list_rlk = re.findall(r"'?\w[\w']*(?:-\w+)*'?'", total_rlk_join) #regex get words all characters
    total_list_rlk_mod = [item.replace("'", "") for item in total_list_rlk]

    total_list_rlk_unique = list(set(total_list_rlk_mod))#items unique
    #total_list_words = '\n'.join(total_list_unique)
    
    with open (filename + "_rlk_total_domains.csv", "w") as f:
        wr = csv.writer(f, delimiter = "\n")
        wr.writerow(["protein_id"])
        wr.writerow(total_list_rlk_unique)
        f.close()
                        
    # create a summary file with results per class     
    with open (filename + "_rlk_general_results.csv", "w") as f:
        rlk_lrr = set_rlk_lrr.count()
        f.write( 'rlk_lrr,' + str(rlk_lrr) + '\n')
        
        rlk_llectin = set_rlk_llectin.count()
        f.write( 'rlk_l-lectine,' + str(rlk_llectin) + '\n')
        
        rlk_clectin = set_rlk_clectin.count()
        f.write( 'rlk_c-lectin,' + str(rlk_clectin) + '\n')
        
        rlk_glectin_b = set_rlk_glectin_b.count()
        f.write( 'rlk_g-lectin-b,' + str(rlp_glectin_b) + '\n')
        
        rlk_glectin_s = set_rlk_glectin_s.count()
        f.write( 'rlk_g-lectin-s,' + str(rlp_glectin_s) + '\n')
        
        rlk_glectin_p = set_rlk_glectin_p.count()
        f.write( 'rlk_g-lectin-p,' + str(rlp_glectin_p) + '\n')
        
        rlk_lysm = set_rlk_lysm.count()
        f.write('rlk_lysm,' + str(rlk_lysm) + '\n')
        
        rlk_malectin = set_rlk_malectin.count()
        f.write('rlk_malectin,' + str(rlk_malectin) + '\n')
        
        rlk_pr5k = set_rlk_pr5k.count()
        f.write('rlk_pr5k,' + str(rlk_pr5k) + '\n')
        
        rlk_tnfr = set_rlk_tnfr.count()
        f.write('rlk_tnfr,' + str(rlk_tnfr) + '\n')
        
        rlk_wak = set_rlk_wak.count()
        f.write('rlk_wak,' + str(rlk_wak) + '\n')
                
        rlk_egf = set_rlk_egf.count()
        f.write('rlk_egf,' + str(rlk_egf) + '\n')
        
        rlk_stress_antifug = set_rlk_str_anti.count()
        f.write('rlk_stress_antifug,' + str(rlk_stress_antifug) + '\n')
        
        rlk_extra = set_rlk_extra_or_rlck.count()
        f.write('rlk_extra_or_rlck,' + str(rlk_extra) + '\n')
        
        #rlck = set_rlck.count()
        #f.write('rlck,' + str(rlck) + '\n')
        
        rlk_or_rlck_only_pkinase = set_rlk_or_rlck_only_pkinase.count()
        f.write('rlk_or_rlck_only_pkinase,' + str(rlk_or_rlck_only_pkinase) + '\n')
        
        rlk_total = (rlk_lrr + rlk_llectin + rlk_clectin + rlk_glectin_b 
                     + rlk_glectin_s + rlk_glectin_p + rlk_lysm + rlk_pr5k 
                     + rlk_wak + rlk_malectin + rlk_egf + rlk_stress_antifug 
                     + rlk_extra + rlk_or_rlck_only_pkinase)
                     
        f.write('rlk_total,' + str(rlk_total) + '\n')
