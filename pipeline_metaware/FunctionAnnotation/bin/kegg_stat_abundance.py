#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-02 14:45:20
@modify: 2022-06-06
'''


## kegg各个层级功能丰度统计 
## 包含nuohe定义的绝对丰度,相对丰度,分组丰度等信息 

import pandas as pd
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

    anno_df = pd.read_csv(args['kegg_anno'], sep='\t')
    abundance_df = pd.read_csv(args['abundance_table'], sep='\t').rename(columns={'gene': 'Query'})
    merge_df = pd.merge(anno_df, abundance_df, on='Query', how='outer').fillna('Others')
    
    ################# 处理ko
    ko_df = merge_df.loc[:, ['Query', 'KO'] + samples]
    ko_df.drop_duplicates(subset=['Query', 'KO'], inplace=True)
    stat_ko = ko_df.groupby('KO').sum()
    stat_ko.loc['Others', :] = stat_ko.loc['Others', :] + stat_ko.loc['-', :]
    stat_ko.drop(index=['-'], inplace=True)
    # 列归一化 * 放大系数
    stat_ko_re = (stat_ko / stat_ko.sum(axis=0))
    stat_ko_ab = (stat_ko / stat_ko.sum(axis=0)) * max(stat_ko.sum(axis=0))
    stat_ko_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.ko.xls'.format(**args), sep='\t', index=None)
    stat_ko_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.ko.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_ko_re_group = get_target_group_df(stat_ko_re, group_info)
    stat_ko_re_group = stat_ko_re_group / stat_ko_re_group.sum(axis=0)
    stat_ko_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.ko.xls'.format(**args), sep='\t', index=None)


    ################# 处理ec
    ec_df = merge_df.loc[:, ['Query', 'KO_EC'] + samples]
    ec_split = ec_df.drop('KO_EC', axis=1).join(ec_df['KO_EC'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('KO_EC'))
    ec_split.drop_duplicates(subset=['Query', 'KO_EC'], inplace=True)  # gene-muti_ko-same_ec
    stat_ec = ec_split.groupby('KO_EC').sum()
    stat_ec.loc['Others', :] = stat_ec.loc['Others', :] + stat_ec.loc['-', :]
    stat_ec.drop(index=['-'], inplace=True)
    # 列归一化 * 放大系数
    stat_ec_re = (stat_ec / stat_ec.sum(axis=0))
    stat_ec_ab = (stat_ec / stat_ec.sum(axis=0)) * max(stat_ec.sum(axis=0))
    stat_ec_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.ec.xls'.format(**args), sep='\t', index=None)
    stat_ec_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.ec.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_ec_re_group = get_target_group_df(stat_ec_re, group_info)
    stat_ec_re_group = stat_ec_re_group / stat_ec_re_group.sum(axis=0)
    stat_ec_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.ec.xls'.format(**args), sep='\t', index=None)


    ################### 处理mo
    mo_df = merge_df.loc[:, ['Query', 'Module'] + samples]
    mo_split = mo_df.drop('Module', axis=1).join(mo_df['Module'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('Module'))
    mo_split.drop_duplicates(subset=['Query', 'Module'], inplace=True)  # gene-muti_ko-same_ec
    stat_mo = mo_split.groupby('Module').sum()
    stat_mo.loc['Others', :] = stat_mo.loc['Others', :] + stat_mo.loc['-', :]
    stat_mo.drop(index=['-'], inplace=True)
    #列归一化 * 放大系数
    stat_mo_re = (stat_mo / stat_mo.sum(axis=0))  
    stat_mo_ab = (stat_mo / stat_mo.sum(axis=0)) * max(stat_mo.sum(axis=0))
    stat_mo_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.module.xls'.format(**args), sep='\t', index=None)
    stat_mo_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.module.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_mo_re_group = get_target_group_df(stat_mo_re, group_info)
    stat_mo_re_group = stat_mo_re_group / stat_mo_re_group.sum(axis=0)
    stat_mo_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.module.xls'.format(**args), sep='\t', index=None)


    #################### 处理pathway
    pathway_df = merge_df.loc[:, ['Query', 'KO_Pathway'] + samples]
    pathway_split = pathway_df.drop('KO_Pathway', axis=1).join(pathway_df['KO_Pathway'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('KO_Pathway'))  
    level_df = pathway_split['KO_Pathway'].str.split(';', expand=True).fillna('-')
    level_df.columns = ['map', 'KO_Pathway_Level1', 'KO_Pathway_Level2', 'KO_Pathway_Level3']
    pathway_list = pd.concat([pathway_split, level_df], axis=1)
    pathway_list['KO_Pathway_Level2'] = pathway_list['KO_Pathway_Level1'] + ';' +pathway_list['KO_Pathway_Level2']  # 改写名称

    ## 层级1
    level1_df = pathway_list.loc[:, ['Query', 'KO_Pathway_Level1'] + samples]
    level1_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level1'], inplace=True)
    stat_level1 = level1_df.groupby('KO_Pathway_Level1').sum()
    stat_level1.loc['Others', :] = stat_level1.loc['-', :]
    stat_level1.drop(index=['-'], inplace=True)
    #列归一化 * 放大系数
    stat_level1_re = (stat_level1 / stat_level1.sum(axis=0)) 
    stat_level1_ab = (stat_level1 / stat_level1.sum(axis=0)) * max(stat_level1.sum(axis=0))
    stat_level1_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.level1.xls'.format(**args), sep='\t', index=None)
    stat_level1_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.level1.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_level1_re_group = get_target_group_df(stat_level1_re, group_info)
    stat_level1_re_group = stat_level1_re_group / stat_level1_re_group.sum(axis=0)
    stat_level1_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.level1.xls'.format(**args), sep='\t', index=None)


    ## 层级2
    level2_df = pathway_list.loc[:, ['Query', 'KO_Pathway_Level2'] + samples]
    level2_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level2'], inplace=True)
    stat_level2 = level2_df.groupby('KO_Pathway_Level2').sum()
    stat_level2.loc['Others', :] = stat_level2.loc['-;-', :]
    stat_level2.drop(index=['-;-'], inplace=True)
    #列归一化 * 放大系数
    stat_level2_re = (stat_level2 / stat_level2.sum(axis=0)) 
    stat_level2_ab = (stat_level2 / stat_level2.sum(axis=0)) * max(stat_level2.sum(axis=0))
    stat_level2_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.level2.xls'.format(**args), sep='\t', index=None)
    stat_level2_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.level2.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_level2_re_group = get_target_group_df(stat_level2_re, group_info)
    stat_level2_re_group = stat_level2_re_group / stat_level2_re_group.sum(axis=0)
    stat_level2_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.level2.xls'.format(**args), sep='\t', index=None)


    ## 层级3
    level3_df = pathway_list.loc[:, ['Query', 'map'] + samples].rename(columns={'map': 'KO_Pathway_Level3'})
    level3_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level3'], inplace=True)
    stat_level3 = level3_df.groupby('KO_Pathway_Level3').sum()
    stat_level3.loc['Others', :] = stat_level3.loc['-', :] + stat_level3.loc['Others', :]
    stat_level3.drop(index=['-'], inplace=True)
    #列归一化 * 放大系数
    stat_level3_re = (stat_level3 / stat_level3.sum(axis=0)) 
    stat_level3_ab = (stat_level3 / stat_level3.sum(axis=0)) * max(stat_level3.sum(axis=0))
    stat_level3_re.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.level3.xls'.format(**args), sep='\t', index=None)
    stat_level3_ab.reset_index().to_csv('{result_dir}/Absolute/{result_suffix}.absolute.level3.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_level3_re_group = get_target_group_df(stat_level3_re, group_info)
    stat_level3_re_group = stat_level3_re_group / stat_level3_re_group.sum(axis=0)
    stat_level3_re_group.reset_index().to_csv('{result_dir}/Relative/{result_suffix}.relative.group.level3.xls'.format(**args), sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='kegg 丰度统计分析')
    parser.add_argument('--kegg_anno', help='基因-kegg功能注释结果(Unigenes.blast.m8.filter.anno.xls)')
    parser.add_argument('--abundance_table', help='物种丰度注释结果(Unigenes.absolute.total.tax.xls)')
    parser.add_argument('--sample_file', help='样本配置文件')
    parser.add_argument('--result_suffix', help='输出文件前缀')
    parser.add_argument('--result_dir', help='输出结果目录')

    args = vars(parser.parse_args())

    main(args)
