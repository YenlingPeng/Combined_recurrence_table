# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 17:08:21 2020

@author: Jacob-Lab
"""
import pandas as pd
import os as os

input_path='D:\\Jacob-Lab\\github\\Combined_recurrence_table\\inputs'
output_path='D:\\Jacob-Lab\\github\\Combined_recurrence_table\\outputs'

snv_recurrence = pd.read_excel("D:\\Jacob-Lab\\github\\SNV_recurrence_analysis_tool\\outputs\\germline_snv_VaraintBasedArray.xlsx", index_col="Idx")
sv_recurrence = pd.read_excel("D:\\Jacob-Lab\\github\\StructureVariant_recurrence_analysis\\outputs\\germline_sv_VaraintBasedArray.xlsx", index_col="Idx")

# change the directory
os.chdir('{}'.format(input_path))
  
# add snv ACMG_ranking info (from TaiGenomics)
ranking_info = pd.read_table("union_of_ranking_from_TaiGenomics.txt")
df1 = snv_recurrence[["Chr","End"]]
df2 = ranking_info[["Chr","End","ACMG_Classes"]].drop_duplicates()
mergedStuff1 = pd.merge(df1, df2, on=['Chr','End'], how='left')
ACMG_Classes = mergedStuff1["ACMG_Classes"].fillna("NaN").tolist()

# Unified field format for recurrence
snv = {'Chr':snv_recurrence["Chr"],
       'Start':snv_recurrence["Start"],
       'End':snv_recurrence["End"],
       'SV_length':"NaN",
       'SV_type':"NaN",
       'Func.refGene/Location1':snv_recurrence["Func.refGene"],
       'ExonicFunc.refGene/location2':snv_recurrence["ExonicFunc.refGene"],
       'Gene_name':snv_recurrence["Gene.refGene"],
       'AAChange':snv_recurrence["AAChange.refGene"],
       'Variant_ID':snv_recurrence["avsnp150"],
       'GD_AF_EAS':snv_recurrence["AF_eas.1"],
       'GD_POPMAX_AF':"NaN",
       'ACMG_classes/AnnotSV_ranking': ACMG_Classes,
       'Variant_type':"snv",
       'germline_test01':snv_recurrence["germline_test01"],
       'germline_test02':snv_recurrence["germline_test02"],
       'Counts':snv_recurrence["Counts"]
       }

snv = pd.DataFrame(snv)

sv = {'Chr':sv_recurrence["SV chrom"],
       'Start':sv_recurrence["SV start"],
       'End':sv_recurrence["SV end"],
       'SV_length':sv_recurrence["SV length"],
       'SV_type':sv_recurrence["SV type"],
       'Func.refGene/Location1':sv_recurrence["location"],
       'ExonicFunc.refGene/location2':sv_recurrence["location2"],
       'Gene_name':sv_recurrence["Gene name"],
       'AAChange':"NaN",
       'Variant_ID':sv_recurrence["GD_ID"],
       'GD_AF_EAS':sv_recurrence["GD_AF"],
       'GD_POPMAX_AF':sv_recurrence["GD_POPMAX_AF"],
       'ACMG_classes/AnnotSV_ranking':sv_recurrence["AnnotSV ranking"],
       'Variant_type':sv_recurrence["Variant_type"],
       'germline_test01':sv_recurrence["germline_test01"],
       'germline_test02':sv_recurrence["germline_test02"],
       'Counts':sv_recurrence["Counts"]
       }
sv = pd.DataFrame(sv)

# combined the snv/sv variant based array
recurrence_table = pd.concat ([snv,sv] , axis = 0)

# candidate gene list info
# "+" = candidate gene & approved 
# "*" = candidate gene & non-approved
# "-" = not candidate gene & approved
# "other" = not candidate gene not approved
# "x" = the gene not in the check list

candidate = pd.read_table("deafness_hgnc_result_20200921.txt")
df3 = pd.DataFrame(recurrence_table["Gene_name"])
df4 = pd.DataFrame(candidate[["Gene_name","filter"]].drop_duplicates())
mergedStuff2 = pd.merge(df3, df4, on=['Gene_name'], how='left').fillna("NaN")
gene_filter =  mergedStuff2["filter"].tolist()

# add candidate gene filter to recurrence table
recurrence_table.insert(14,"Candidate_gene_filter" , gene_filter) 
recurrence_table.to_excel('{}'.format(output_path) + "\\germline_combined_recurrence_table.xlsx", index= False)

