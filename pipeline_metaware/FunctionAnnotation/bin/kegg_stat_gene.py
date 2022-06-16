#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-01 17:13:09
'''


# 统计不同KEGG 不同水平下基因数目统计 

import pandas as pd  
from collections import defaultdict


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



def main(args):  # sourcery skip: low-code-quality
    '''
    {
        ko: {
            'ko1': [],  
            'ko2': []
        },
        
        'ec': {

        }
    }
    '''
    stat_info = defaultdict(dict)
    classes = ['ko', 'ec', 'module', 'level1', 'level2', 'level3']
    for _class in classes:
        stat_info[_class] = defaultdict(list)

    anno_df = pd.read_csv(args['kegg_anno'], sep='\t')
    gene_df = pd.read_csv(args['gene_table'], sep='\t').rename(columns={'gene': 'Query'})
    merge_df = pd.merge(anno_df, gene_df, on='Query', how='outer').fillna('Others')

    headers = list(merge_df.columns)
    header_info = get_index(headers)

    for row in merge_df.values:
        query = row[header_info['query']] 
        ko = row[header_info['ko']]
        ec = row[header_info['ko_ec']]
        pathway = row[header_info['ko_pathway']]
        mo = row[header_info['module']]

        # 处理ko, 一个基因唯一ko
        if ko.startswith(('-', 'Others')):
            stat_info['ko']['Others'].append(query)
        else:
            stat_info['ko'][ko].append(query)  

        # 处理ec, 一个ko,可能对应多个ec
        if ec.startswith(('-', 'Others')):
            stat_info['ec']['Others'].append(query)
        else:
            eclist = ec.split('|')
            for _ec in eclist:
                stat_info['ec'][_ec].append(query)

        # 处理mo
        if mo.startswith(('-', 'Others')):
            stat_info['module']['Others'].append(query)
        else:
            molist = mo.split('|')
            for _mo in molist:
                stat_info['module'][_mo].append(query)

        # 处理pathway,除了组合,还需要分成3个level,
        # 拆分过程中,不同的map可能具有相同的level和level2
        if pathway.startswith(('-', 'Others')):
            stat_info['level1']['Others'].append(query)
            stat_info['level2']['Others'].append(query)
            stat_info['level3']['Others'].append(query)
        else:
            pathway_list = pathway.split('|')
            for p in pathway_list:
                pp = p.split(';')
                _3 = pp[0]
                _1 = pp[1]
                _2 = ';'.join(pp[1:3])

                stat_info['level1'][_1].append(query)
                stat_info['level2'][_2].append(query)
                stat_info['level3'][_3].append(query)


    for _class in stat_info:
        out = '{result_suffix}.combine.{_class}.xls'.format(**args, **locals())
        with open(out, 'w') as fw:
            for item in stat_info[_class]:  # item: K0001
                ## 需要去重,
                ## 1. 一个KO多个map时,不同map在level1和level2上可能一样
                ## 2. 一个基因映射多个KO, 这多个KO可能映射相同的ec或者module。
                gene_num = len(set(stat_info[_class][item]))  
                gene_ids = ','.join(set(stat_info[_class][item]))
                fw.write('{item}\t{gene_num}\t{gene_ids}\n'.format(**locals()))


    #################################################################################
    ###################【统计每个层级，每个样本的基因数】##############################
    samples = get_sampels(args['sample_file'])
    for sample in samples:
        merge_df.loc[merge_df.loc[:, sample] > 0, sample] = 1  # 转换为1,方便统计基因数
    
    ## 处理ko
    ko_df = merge_df.loc[:, ['Query', 'KO'] + samples]
    ko_df.drop_duplicates(subset=['Query', 'KO'], inplace=True)
    stat_ko = ko_df.groupby('KO').sum()
    stat_ko.loc['Others', :] = stat_ko.loc['Others', :] + stat_ko.loc['-', :]
    stat_ko.drop(index=['-'], inplace=True)
    stat_ko.reset_index().to_csv('{result_suffix}.ko.xls'.format(**args), sep='\t', index=None)

    ## 处理ec
    ec_df = merge_df.loc[:, ['Query', 'KO_EC'] + samples]
    ec_split = ec_df.drop('KO_EC', axis=1).join(ec_df['KO_EC'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('KO_EC'))
    ec_split.drop_duplicates(subset=['Query', 'KO_EC'], inplace=True)  # gene-muti_ko-same_ec
    stat_ec = ec_split.groupby('KO_EC').sum()
    stat_ec.loc['Others', :] = stat_ec.loc['Others', :] + stat_ec.loc['-', :]
    stat_ec.drop(index=['-'], inplace=True)
    stat_ec.reset_index().to_csv('{result_suffix}.ec.xls'.format(**args), sep='\t', index=None)

    ## 处理mo
    mo_df = merge_df.loc[:, ['Query', 'Module'] + samples]
    mo_split = mo_df.drop('Module', axis=1).join(mo_df['Module'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('Module'))
    mo_split.drop_duplicates(subset=['Query', 'Module'], inplace=True)  # gene-muti_ko-same_ec
    stat_mo = mo_split.groupby('Module').sum()
    stat_mo.loc['Others', :] = stat_mo.loc['Others', :] + stat_mo.loc['-', :]
    stat_mo.drop(index=['-'], inplace=True)
    stat_mo.reset_index().to_csv('{result_suffix}.module.xls'.format(**args), sep='\t', index=None)
    
    ## 处理pathway
    pathway_df = merge_df.loc[:, ['Query', 'KO_Pathway'] + samples]
    pathway_split = pathway_df.drop('KO_Pathway', axis=1).join(pathway_df['KO_Pathway'].str.split('|', expand=True).stack().reset_index(level=1, drop=True).rename('KO_Pathway'))  
    level_df = pathway_split['KO_Pathway'].str.split(';', expand=True).fillna('-')
    level_df.columns = ['map', 'KO_Pathway_Level1', 'KO_Pathway_Level2', 'KO_Pathway_Level3']
    pathway_list = pd.concat([pathway_split, level_df], axis=1)
    pathway_list['KO_Pathway_Level2'] = pathway_list['KO_Pathway_Level1'] + ';' +pathway_list['KO_Pathway_Level2']  # 改写名称

    # 层级1
    level1_df = pathway_list.loc[:, ['Query', 'KO_Pathway_Level1'] + samples]
    level1_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level1'], inplace=True)
    stat_level1 = level1_df.groupby('KO_Pathway_Level1').sum()
    stat_level1.loc['Others', :] = stat_level1.loc['-', :]
    stat_level1.drop(index=['-'], inplace=True)
    stat_level1.reset_index().to_csv('{result_suffix}.level1.xls'.format(**args), sep='\t', index=None)

    # 层级2
    level2_df = pathway_list.loc[:, ['Query', 'KO_Pathway_Level2'] + samples]
    level2_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level2'], inplace=True)
    stat_level2 = level2_df.groupby('KO_Pathway_Level2').sum()
    stat_level2.loc['Others;Others', :] = stat_level2.loc['-;-', :]
    stat_level2.drop(index=['-;-'], inplace=True)
    stat_level2.reset_index().to_csv('{result_suffix}.level2.xls'.format(**args), sep='\t', index=None)


    # 层级3
    level3_df = pathway_list.loc[:, ['Query', 'map'] + samples].rename(columns={'map': 'KO_Pathway_Level3'})
    level3_df.drop_duplicates(subset=['Query', 'KO_Pathway_Level3'], inplace=True)
    stat_level3 = level3_df.groupby('KO_Pathway_Level3').sum()
    stat_level3.loc['Others', :] = stat_level3.loc['-', :] + stat_level3.loc['Others', :]
    stat_level3.drop(index=['-'], inplace=True)
    stat_level3.reset_index().to_csv('{result_suffix}.level3.xls'.format(**args), sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='基因数目统计')
    parser.add_argument('--kegg_anno', help='kegg 注释结果')
    parser.add_argument('--gene_table', help='样本基因计数结果')
    parser.add_argument('--sample_file', help='sample file, 提取样本名称')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())

    main(args)

