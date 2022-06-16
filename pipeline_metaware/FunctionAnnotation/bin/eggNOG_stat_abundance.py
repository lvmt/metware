#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-15 15:34:51
'''

# same method just like before 
from random import sample
import pandas as pd   
from collections import defaultdict
import re  


_class_info = {
    'J':  'J: Translation, ribosomal structure and biogenesis',
    'A':  'A: RNA processing and modification',
    'K':  'K: Transcription ',
    'L':  'L: Replication, recombination and repair' ,
    'B':  'B: Chromatin structure and dynamics' ,
    'D':  'D: Cell cycle control, cell division, chromosome partitioning',
    'Y':  'Y: Nuclear structure',
    'V':  'V: Defense mechanisms' ,
    'T':  'T: Signal transduction mechanisms' ,
    'M':  'M: Cell wall/membrane/envelope biogenesis' ,
    'N':  'N: Cell motility' ,
    'Z':  'Z: Cytoskeleton ',
    'W':  'W: Extracellular structures' ,
    'U':  'U: Intracellular trafficking, secretion, and vesicular transport ',
    'O':  'O: Posttranslational modification, protein turnover, chaperones' ,
    'C':  'C: Energy production and conversion' ,
    'G':  'G: Carbohydrate transport and metabolism' ,
    'E':  'E: Amino acid transport and metabolism' ,
    'F':  'F: Nucleotide transport and metabolism ',
    'H':  'H: Coenzyme transport and metabolism' ,
    'I':  'I: Lipid transport and metabolism' ,
    'P':  'P: Inorganic ion transport and metabolism' ,
    'Q':  'Q: Secondary metabolites biosynthesis, transport and catabolism' ,
    'R':  'R: General function prediction only' ,
    'S':  'S: Function unknown',
    'Others': 'Others'
}


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

    anno_df = pd.read_csv(args['eggNOG_anno'], sep='\t').fillna('Others')
    abundance_df = pd.read_csv(args['abundance_table'], sep='\t').rename(columns={'gene': 'query'})
    merge_df = pd.merge(anno_df, abundance_df, on='query', how='outer').fillna('Others')
     
    ####### 处理_class（level1）
    class_df = merge_df.loc[:, ['query', 'Functional_Category'] + samples]
    class_df['Functional_Category'] = class_df['Functional_Category'].apply(lambda x: '|'.join(list(x)) if not x.startswith('Others') else x ) #有些对应多个_class
    class_split = class_df.drop('Functional_Category', axis=1).join(class_df['Functional_Category'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('Functional_Category'))
    class_split.drop_duplicates(subset=['query', 'Functional_Category'], inplace=True)
    class_split['Functional_Category'] = class_split['Functional_Category'].map(_class_info)
    stat_class = class_split.groupby('Functional_Category').sum()

    stat_class_re = stat_class / stat_class.sum(axis=0)
    stat_class_ab = (stat_class / stat_class.sum(axis=0)) * max(stat_class.sum(axis=0))
    stat_class_re.reset_index().to_csv('{result_suffix}.relative.level1.xls'.format(**args), sep='\t', index=None)
    stat_class_ab.reset_index().to_csv('{result_suffix}.absolute.level1.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_class_group = get_target_group_df(stat_class_re, group_info)
    stat_class_group = stat_class_group / stat_class_group.sum(axis=0)
    stat_class_group.reset_index().to_csv('{result_suffix}.relative.group.level1.xls'.format(**args), sep='\t', index=None)


    #### 处理desc，（level2）
    level2_df = merge_df.loc[:, ['query', 'OG_Description'] + samples]
    level2_df.drop_duplicates(subset=['query', 'OG_Description'], inplace=True)
    stat_level2 = level2_df.groupby('OG_Description').sum()

    stat_level2_re = stat_level2 / stat_level2.sum(axis=0)
    stat_level2_ab = (stat_level2 / stat_level2.sum(axis=0)) * max(stat_level2.sum(axis=0))
    stat_level2_re.reset_index().to_csv('{result_suffix}.relative.level2.xls'.format(**args), sep='\t', index=None)
    stat_level2_ab.reset_index().to_csv('{result_suffix}.absolute.level2.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_level2_group = get_target_group_df(stat_level2_re, group_info)
    stat_level2_group = stat_level2_group / stat_level2_group.sum(axis=0)
    stat_level2_group.reset_index().to_csv('{result_suffix}.relative.group.level2.xls'.format(**args), sep='\t', index=None)


    #### 处理OG
    og_df = merge_df.loc[:, ['query', 'Ortholog_Group'] + samples]
    og_df.drop_duplicates(subset=['query', 'Ortholog_Group'], inplace=True)
    stat_og = og_df.groupby('Ortholog_Group').sum()

    stat_og_re = stat_og / stat_og.sum(axis=0)
    stat_og_ab = (stat_og / stat_og.sum(axis=0)) * max(stat_og.sum(axis=0))
    stat_og_re.reset_index().to_csv('{result_suffix}.relative.og.xls'.format(**args), sep='\t', index=None)
    stat_og_ab.reset_index().to_csv('{result_suffix}.absolute.og.xls'.format(**args), sep='\t', index=None)
    # 分组结果
    stat_og_group = get_target_group_df(stat_og_re, group_info)
    stat_og_group = stat_og_group / stat_og_group.sum(axis=0)
    stat_og_group.reset_index().to_csv('{result_suffix}.relative.group.og.xls'.format(**args), sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='kegg 丰度统计分析')
    parser.add_argument('--eggNOG_anno', help='基因-eggNOG功能注释结果(Unigenes.blast.m8.filter.anno.xls)')
    parser.add_argument('--abundance_table', help='物种丰度注释结果(Unigenes.absolute.total.tax.xls)')
    parser.add_argument('--sample_file', help='样本配置文件')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())

    main(args)
