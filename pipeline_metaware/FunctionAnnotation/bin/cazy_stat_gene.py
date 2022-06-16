#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-14 11:40:55
'''


# 统计不同CAZy 不同水平下基因数目统计 
import pandas as pd  
from collections import defaultdict
import re


def get_index(headers):
    if not isinstance(headers, list):
        headers = headers.strip('\n').split('\t')

    return {item.lower(): index for index,item in enumerate(headers)}


def get_sampels(sample_file):
    samples = []
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            sample = linelist[0]
            samples.append(sample)
    return samples


def main(args):
    anno_df = pd.read_csv(args['cazy_anno'], sep='\t').fillna('Others')
    gene_df = pd.read_csv(args['gene_table'], sep='\t').rename(columns={'gene': 'query'})
    merge_df = pd.merge(anno_df, gene_df, on='query', how='outer').fillna('Others')

    headers = list(merge_df.columns)
    header_info = get_index(headers)

    
    samples = get_sampels(args['sample_file'])
    for sample in samples:
        merge_df.loc[merge_df.loc[:, sample] > 0, sample] = 1  # 转换为1,方便统计基因数

    ##################### 处理6大通路信息 #############################
    cazy_df = merge_df.loc[:, ['query', 'family'] + samples]
    cazy_df.rename(columns={'family': 'CAZy_Family'}, inplace=True)
    cazy_split = cazy_df.drop('CAZy_Family', axis=1).join(cazy_df['CAZy_Family'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('CAZy_Family'))
    cazy_split['CAZy_Modules'] = cazy_split['CAZy_Family'].apply(lambda x: re.search('[a-zA-Z]+', x).group())
    
    ## 处理level1
    cazy_level1_df = cazy_split.loc[:, ['query', 'CAZy_Modules'] + samples]
    cazy_level1_df.drop_duplicates(subset=['query', 'CAZy_Modules'], inplace=True)
    stat_cazy_level1 = cazy_level1_df.groupby('CAZy_Modules').sum()
    stat_cazy_level1.reset_index().to_csv('{result_suffix}.level1.xls'.format(**args), sep='\t', index=None)

    # 汇总统计
    level1_info = defaultdict(set)
    for row in cazy_level1_df.loc[:, ['query', 'CAZy_Modules']].values:
        gene, _class = row 
        level1_info[_class].add(gene)

    with open('{result_suffix}.combine.level1.xls'.format(**args), 'w') as fw:
        fw.write('CAZy_Modules\tGene_Num\tGene_IDs\n')
        for _class in level1_info:
            gene_nums = len(level1_info[_class])
            gene_ids = ','.join(level1_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')


    ## 处理level2
    cazy_level2_df = cazy_split.loc[:, ['query', 'CAZy_Family'] + samples]
    cazy_level2_df['CAZy_Family'] = cazy_level2_df['CAZy_Family'].apply(lambda x: re.search(r'[a-zA-Z0-9]+', x).group())  # 将亚家族转换为家族
    cazy_level2_df.drop_duplicates(subset=['query', 'CAZy_Family'], inplace=True)
    stat_cazy_level2 = cazy_level2_df.groupby('CAZy_Family').sum()
    stat_cazy_level2.reset_index().to_csv('{result_suffix}.level2.xls'.format(**args), sep='\t', index=None)

    # 汇总统计
    level2_info = defaultdict(set)
    for row in cazy_level2_df.loc[:, ['query', 'CAZy_Family']].values:
        gene, _class = row 
        level2_info[_class].add(gene)

    with open('{result_suffix}.combine.level2.xls'.format(**args), 'w') as fw:
        fw.write('CAZy_Family\tGene_Num\tGene_IDs\n')
        for _class in level2_info:
            gene_nums = len(level2_info[_class])
            gene_ids = ','.join(level2_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')


    ############################ 处理EC 
    ec_df = merge_df.loc[:, ['query', 'ec'] + samples]
    ec_split = ec_df.drop('ec', axis=1).join(ec_df['ec'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('ec'))
    ec_split.drop_duplicates(subset=['query', 'ec'], inplace=True)
    stat_ec = ec_split.groupby('ec').sum()
    stat_ec.reset_index().to_csv('{result_suffix}.ec.xls'.format(**args), sep='\t', index=None)

    # 汇总统计
    ec_info = defaultdict(set)
    for row in ec_split.loc[:, ['query', 'ec']].values:
        gene, _class = row 
        ec_info[_class].add(gene)

    with open('{result_suffix}.combine.ec.xls'.format(**args), 'w') as fw:
        fw.write('EC\tGene_Num\tGene_IDs\n')
        for _class in ec_info:
            gene_nums = len(ec_info[_class])
            gene_ids = ','.join(ec_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='基因数目统计')
    parser.add_argument('--cazy_anno', help='cazy 注释结果')
    parser.add_argument('--gene_table', help='样本基因计数结果')
    parser.add_argument('--sample_file', help='sample file, 提取样本名称')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())

    main(args)

