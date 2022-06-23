#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-06 09:59:46
'''


# 合并全部样本的丰度结果
# 合并全部样本的基因数目结果(后续相对转换为绝对定量的时候需要)
import pandas as pd  
import os  
from functools import reduce



import argparse
parser = argparse.ArgumentParser(description='merge abundance table')
parser.add_argument('--abundance_tables', nargs='+', help='丰度文件')
parser.add_argument('--merge_result_suffix', help='输出结果')

args = vars(parser.parse_args())
# print(args)

df_abundance_list = []
df_gene_list = []
for table in args['abundance_tables']:
    name = os.path.basename(table).replace('.abundance.table', '')
    df = pd.read_csv(table, sep='\t')
    df_abundance = df.loc[:, ['gene', 'abundance']].rename(columns={'abundance': name})
    df_gene = df.loc[:, ['gene', 'count']].rename(columns={'count': name})
    df_abundance_list.append(df_abundance)
    df_gene_list.append(df_gene)


combine_abundance_df = reduce(lambda x,y: pd.merge(x, y, on='gene', how='outer').fillna(0), df_abundance_list)
combine_gene_df = reduce(lambda x,y: pd.merge(x, y, on='gene', how='outer').fillna(0), df_gene_list)

merge_gene = '{merge_result_suffix}.gene.readsNum'.format(**args)
merge_abudance = '{merge_result_suffix}.abundance.relative.xls'.format(**args)

combine_abundance_df.to_csv(merge_abudance, sep='\t', index=None)
combine_gene_df.to_csv(merge_gene, sep='\t', index=None)

