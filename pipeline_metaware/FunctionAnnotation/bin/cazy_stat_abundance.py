#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-14 09:35:29
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-14 09:35:29
'''

# cazy数据库，统计丰度信息；思路同kegg


from tkinter.messagebox import NO
import pandas as pd
import re 
from collections import defaultdict


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


def get_group_info(sample_file):
    # 获取样本的分组信息
    group_info = defaultdict(list)
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            sample = linelist[0]
            group = linelist[1]
            group_info[group].append(sample)
    return group_info  


def get_target_group_df(df, group_info):
    # 获取分组的dataframe
    for _group in group_info:
        df[_group] = df.loc[:, group_info[_group]].sum(axis=1)

    df = df.loc[:, list(group_info.keys())]
    return df 


def main(args):
    samples = get_sampels(args['sample_file'])
    group_info = get_group_info(args['sample_file'])


    anno_df = pd.read_csv(args['cazy_anno'], sep='\t').fillna('Others')
    abundance_df = pd.read_csv(args['abundance_table'], sep='\t').rename(columns={'gene': 'query'})
    merge_df = pd.merge(anno_df, abundance_df, on='query', how='outer').fillna('Others')


    ########################### 处理6大通路文件
    cazy_df = merge_df.loc[:, ['query', 'family'] + samples]
    cazy_df.rename(columns={'family': 'CAZy_Family'}, inplace=True)
    cazy_split = cazy_df.drop('CAZy_Family', axis=1).join(cazy_df['CAZy_Family'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('CAZy_Family'))
    cazy_split['CAZy_Modules'] = cazy_split['CAZy_Family'].apply(lambda x: re.search('[a-zA-Z]+', x).group())

    ## level1水平
    cazy_level1_df = cazy_split.loc[:, ['query', 'CAZy_Modules'] + samples]
    cazy_level1_df.drop_duplicates(subset=['query', 'CAZy_Modules'], inplace=True)
    stat_cazy_level1 = cazy_level1_df.groupby('CAZy_Modules').sum()
    
    stat_cazy_level1_re = stat_cazy_level1 / stat_cazy_level1.sum(axis=0)
    stat_cazy_level1_ab = (stat_cazy_level1 / stat_cazy_level1.sum(axis=0)) * max(stat_cazy_level1.sum(axis=0))
    stat_cazy_level1_re.reset_index().to_csv('{result_suffix}.relative.level1.xls'.format(**args), sep='\t', index=None)
    stat_cazy_level1_ab.reset_index().to_csv('{result_suffix}.absolute.level1.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_cazy_level1_group = get_target_group_df(stat_cazy_level1_re, group_info)
    stat_cazy_level1_group = stat_cazy_level1_group / stat_cazy_level1_group.sum(axis=0)
    stat_cazy_level1_group.reset_index().to_csv('{result_suffix}.relative.group.level1.xls'.format(**args), sep='\t', index=None)

    ## level2 水平
    cazy_level2_df = cazy_split.loc[:, ['query', 'CAZy_Family'] + samples]
    cazy_level2_df['CAZy_Family'] = cazy_level2_df['CAZy_Family'].apply(lambda x: re.search(r'[a-zA-Z0-9]+', x).group())  # 将亚家族转换为家族
    cazy_level2_df.drop_duplicates(subset=['query', 'CAZy_Family'], inplace=True)
    stat_cazy_level2 = cazy_level2_df.groupby('CAZy_Family').sum()
    
    stat_cazy_level2_re = stat_cazy_level2 / stat_cazy_level2.sum(axis=0)
    stat_cazy_level2_ab = (stat_cazy_level2 / stat_cazy_level2.sum(axis=0)) * max(stat_cazy_level2.sum(axis=0))
    stat_cazy_level2_re.reset_index().to_csv('{result_suffix}.relative.level2.xls'.format(**args), sep='\t', index=None)
    stat_cazy_level2_ab.reset_index().to_csv('{result_suffix}.absolute.level2.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_cazy_level2_group = get_target_group_df(stat_cazy_level2_re, group_info)
    stat_cazy_level2_group = stat_cazy_level2_group / stat_cazy_level2_group.sum(axis=0)
    stat_cazy_level2_group.reset_index().to_csv('{result_suffix}.relative.group.level2.xls'.format(**args), sep='\t', index=None)


    # ################### 处理酶相关
    ec_df = merge_df.loc[:, ['query', 'ec'] + samples]
    ec_split = ec_df.drop('ec', axis=1).join(ec_df['ec'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('ec'))

    ec_split.drop_duplicates(subset=['query', 'ec'], inplace=True)
    ec_stat = ec_split.groupby('ec').sum()
    
    ec_stat_re = ec_stat / ec_stat.sum(axis=0)
    ec_stat_ab = (ec_stat / ec_stat.sum(axis=0)) * max(ec_stat.sum(axis=0))
    ec_stat_re.reset_index().to_csv('{result_suffix}.relative.ec.xls'.format(**args), sep='\t', index=None)
    ec_stat_ab.reset_index().to_csv('{result_suffix}.absolute.ec.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    ec_stat_group = get_target_group_df(ec_stat_re, group_info)
    ec_stat_group = ec_stat_group / ec_stat_group.sum(axis=0)
    ec_stat_group.reset_index().to_csv('{result_suffix}.relative.group.ec.xls'.format(**args), sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='kegg 丰度统计分析')
    parser.add_argument('--cazy_anno', help='基因-kegg功能注释结果(Unigenes.blast.m8.filter.anno.xls)')
    parser.add_argument('--abundance_table', help='物种丰度注释结果(Unigenes.absolute.total.tax.xls)')
    parser.add_argument('--sample_file', help='样本配置文件')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())

    main(args)