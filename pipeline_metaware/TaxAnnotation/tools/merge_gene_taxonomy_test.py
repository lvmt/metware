#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-25 16:33:47
code check: 2022-06-22
'''

'''
目的：
汇总 基因 与 物种注释 
某个物种在某个样品中的基因数目，等于注释为该物种的基因中，丰度不为 0 的基因数目加和
'''


import pandas as pd  
from collections import defaultdict


def get_samples(sample_file):
    samples = []
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            samples.append(linelist[0])
    return samples  


def main(args):
    samples = get_samples(args['sample_file'])
    tax_df = pd.read_csv(args['taxonomy_table'], sep='\t', header=None).rename(columns={0: 'Reference_ID', 1: 'taxonomy'})
    gene = pd.read_csv(args['gene_table'], sep='\t')
    merge_df = pd.merge(tax_df, gene, on='Reference_ID', how='outer').fillna('k__Others;p__Others;c__Others;o__Others;f__Others;g__Others;s__Others')
    merge_df['taxonomy'] = merge_df['taxonomy'].apply(lambda x: x if len(x.split(';')) == 7 else ';'.join(x.split(';') + ['Others'] * (7 - len(x.split(';')))))
    # 因为是统计基因的数目,因此我们需要把样本中大于0的数修改为1,便于后续的分组统计.
    for sample in samples:
        merge_df.loc[merge_df.loc[:, sample] > 0, sample] = 1

    ## 进行层级拆分(每个样本)
    class_info = {
        'kingdom': 1,
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'genus': 6,
        'species': 7
    }

    for _class in class_info:
        merge_df['Taxonomy'] = merge_df['taxonomy'].apply(lambda x:';'.join(x.split(';')[:class_info[_class]]))
        merge_df.loc[merge_df['Taxonomy'].apply(lambda x: 'Unclassified' in x.split(';')[-1]), 'Taxonomy'] = 'Others'
        merge_df.loc[merge_df['Taxonomy'].apply(lambda x: 'Other' in x.split(';')[-1]), 'Taxonomy'] = 'Others'
        stat_df = merge_df.groupby('Taxonomy').sum().reset_index()

        # stat_df[_class] = stat_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[-2:])) 
        stat_df['Taxonomy'] = stat_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[-2:])) 
        out = '{result_suffix}.{_class}.xls'.format(**args, **locals())
        stat_df.to_csv(out, sep='\t', index=None)


        ## 进行层级拆分(全部样本合并统计)
        out = '{result_suffix}.total.{_class}.xls'.format(**args, **locals())
        # merge_df['Detail_Taxonomy'] = merge_df['taxonomy'].apply(lambda x: ';'.join(x.split(';')[:class_info[_class]]))
        _class_info = defaultdict(set)

        for row in merge_df.loc[:, ['Reference_ID', 'Taxonomy']].values:
            gene, tax = row 
            _class_info[tax].add(gene)

        with open(out, 'w') as fw:
            fw.write(f'{_class}\tDetail_Taxonomy\tGene_Num\tGene_IDs\n')
            for _class in _class_info:
                gene_nums = len(_class_info[_class])
                gene_ids = ','.join(_class_info[_class])
                family = '|'.join(_class.split(';')[-2:])
                fw.write(f'{family}\t{_class}\t{gene_nums}\t{gene_ids}\n')

        del _class_info



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='汇总')
    parser.add_argument('--gene_table', help='基因统计表格')
    parser.add_argument('--taxonomy_table', help='物种注释表格')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--result_suffix', help='输出结果前缀')

    args = vars(parser.parse_args())
    main(args)



