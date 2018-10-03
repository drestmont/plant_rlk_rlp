#!/bin/bash
python 9.10_df_pfam_scan_analysis.py &&
python 10.2_df_rlk_rlp_unique_domains.py &&
python 11.1_df_rlk_rlp_parser.py &&
python 12.1_rlk_rlp_merge_files.py &&
python 13.1_rlk_rlp_modified_general_results.py &&
python 14.1_rlk_rlp_report_results.py &&
python 15_rlk_rlp_reindex_order.py
