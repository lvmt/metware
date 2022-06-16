#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-16 10:20:05
'''


import pandas as pd  
from collections import defaultdict
import re


class_info = {
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
    anno_df = pd.read_csv(args['eggNOG_anno'], sep='\t').fillna('Others')
    gene_df = pd.read_csv(args['gene_table'], sep='\t').rename(columns={'Reference_ID': 'Query_id'})
    merge_df = pd.merge(anno_df, gene_df, on='Query_id', how='outer').fillna('Others')

    headers = list(merge_df.columns)
    header_info = get_index(headers)
 
    samples = get_sampels(args['sample_file'])
    for sample in samples:
        merge_df.loc[merge_df.loc[:, sample] > 0, sample] = 1  # 转换为1,方便统计基因数

    #### 处理level1
    level1_df = merge_df.loc[:, ['Query_id', 'Functional_Category'] + samples]
    level1_df['Functional_Category'] = level1_df['Functional_Category'].apply(lambda x: '|'.join(list(x)) if not x.startswith('Others') else x ) #有些对应多个_class
    level1_split = level1_df.drop('Functional_Category', axis=1).join(level1_df['Functional_Category'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('Functional_Category'))
    level1_split.drop_duplicates(subset=['Query_id', 'Functional_Category'], inplace=True)
    level1_split['Functional_Category'] = level1_split['Functional_Category'].map(class_info)

    stat_level1 = level1_split.groupby('Functional_Category').sum()
    stat_level1.reset_index().to_csv('{result_suffix}.level1.xls'.format(**args), sep='\t', index=None)

    ## 汇总统计
    level1_info = defaultdict(set)
    for row in level1_split.loc[:, ['Query_id', 'Functional_Category']].values:
        gene, _class = row
        level1_info[_class].add(gene)

    with open('{result_suffix}.combine.level1.xls'.format(**args), 'w') as fw:
        fw.write('Functional_Category\tGene_Num\tGene_IDs\n')
        for _class in level1_info:
            gene_nums = len(level1_info[_class])
            gene_ids = ','.join(level1_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')


    #### 处理level2
    level2_df = merge_df.loc[:, ['Query_id', 'OG_Description'] + samples]
    level2_df.drop_duplicates(subset=['Query_id', 'OG_Description'], inplace=True)
    
    stat_level2 = level2_df.groupby('OG_Description').sum()
    stat_level2.reset_index().to_csv('{result_suffix}.level2.xls'.format(**args), sep='\t', index=None)

    ## 汇总统计 
    level2_info = defaultdict(set)
    for row in level2_df.loc[:, ['Query_id', 'OG_Description']].values:
        gene, _class = row
        level2_info[_class].add(gene)

    with open('{result_suffix}.combine.level2.xls'.format(**args), 'w') as fw:
        fw.write('OG_Description\tGene_Num\ttGene_IDs\n')
        for _class in level2_info:
            gene_nums = len(level2_info[_class])
            gene_ids = ','.join(level2_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')
            

    #### 处理OG 
    og_df = merge_df.loc[:, ['Query_id', 'Ortholog_Group'] + samples]
    og_df.drop_duplicates(subset=['Query_id', 'Ortholog_Group'], inplace=True)

    stat_og = og_df.groupby('Ortholog_Group').sum()
    stat_og.reset_index().to_csv('{result_suffix}.og.xls'.format(**args), sep='\t', index=None)

    ## 汇总统计
    og_info = defaultdict(set)
    for row in og_df.loc[:, ['Query_id', 'Ortholog_Group']].values:
        gene, _class = row
        og_info[_class].add(gene)

    with open('{result_suffix}.combine.og.xls'.format(**args), 'w') as fw:
        fw.write('Ortholog_Group\tGene_Num\ttGene_IDs\n')
        for _class in og_info:
            gene_nums = len(og_info[_class])
            gene_ids = ','.join(og_info[_class])
            fw.write(f'{_class}\t{gene_nums}\t{gene_ids}\n')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='基因数目统计')
    parser.add_argument('--eggNOG_anno', help='eggNOG 注释结果')
    parser.add_argument('--gene_table', help='样本基因计数结果')
    parser.add_argument('--sample_file', help='sample file, 提取样本名称')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())

    main(args)