#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-06 09:59:46
'''


# 合并全部样本的丰度结果
import pandas as pd  
import os  




import argparse
parser = argparse.ArgumentParser(description='merge abundance table')
parser.add_argument('--abundance_tables', nargs='+', help='丰度文件')
parser.add_argument('--merge_abundance', help='输出结果')

args = vars(parser.parse_args())
# print(args)

df_list = []
for table in args['abundance_tables']:
    name = os.path.basename(table)
    df = pd.read_csv(table, sep='\t')
    df = df.loc[:, ['gene', 'abundance']].rename(columns={'abundance': name})
    df_list.append(df)


combine_df = df_list[0]

for df in df_list[1:]:
    combine_df = pd.merge(combine_df, df, on='gene', how='outer').fillna(0)


combine_df.to_csv(args['merge_abundance'], sep='\t', index=None)

