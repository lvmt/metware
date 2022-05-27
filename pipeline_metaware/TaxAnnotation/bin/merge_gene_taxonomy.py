#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-25 16:33:47
'''



import contextlib
from tokenize import group
# 汇总基因 与 物种注释表格 
# 于某个物种在某个样品中的基因数目，等于在注释为该物种的基因中，丰度不为 0 的基因数目

import pandas as pd  


def main(args):
    tax = pd.read_csv(args['taxonomy_tabel'], sep='\t')
    tax = tax.drop(columns=['taxid'])
    gene = pd.read_csv(args['gene_table'], sep='\t')

    merge_df = pd.merge(tax, gene, on='gene', how='outer').fillna('k__Other|p__Other|c__Other|o__Other|f__Other|g__Other|s__Other')
    # 因为是统计基因的数目,因此我们需要把样本中大于0的数修改为1,便于后续的分组统计.
    samples = list(merge_df.columns)
    with contextlib.suppress(Exception):
        samples.remove('taxonomy')
        samples.remove('gene') 

    for sample in samples:
        merge_df.loc[merge_df.loc[:, sample] > 0, sample] = 1

    # print(merge_df)

    ## 进行层级拆分(每个样本)
    class_info = {
        'kingdom': 1,
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'Genus': 6,
        'species': 7
    }
    for _class in class_info:
        merge_df['Taxonomy'] = merge_df['taxonomy'].apply(lambda x:'|'.join(x.split('|')[:class_info[_class]]))
        class_df = merge_df.groupby('Taxonomy').agg('sum').reset_index()
        class_df[_class] = class_df['Taxonomy'].apply(lambda x:'|'.join(x.split('|')[-2:])) 

        out = '{result_suffix}.{_class}.xls'.format(**args, **locals())
        class_df.to_csv(out, sep='\t', index=None)


    ## 进行层级拆分(全部样本合并)
    for _class in class_info:

        out = '{result_suffix}.combine_samples.{_class}.xls'.format(**args, **locals())
        merge_df['Taxonomy'] = merge_df['taxonomy'].apply(lambda x:'|'.join(x.split('|')[:class_info[_class]]))
        class_df_groups = list(merge_df.groupby('Taxonomy'))

        with open(out, 'w') as fw:
            fw.write(f'{_class}\tdetail_taxonomy\tgene_num\tgene_ids\n')
            for _group in class_df_groups:
                group_name = _group[0]
                simple_group_name = '|'.join(group_name.split('|')[-2:])
                group_df = _group[1]
                group_df['count'] = group_df.sum(axis=1)  # 统计每个微生物在全部样本中出现的次数
                genes = group_df.loc[group_df['count'] > 0, 'gene'].values
                genes_num = len(genes)
                genes = '|'.join(genes)
                fw.write(f'{simple_group_name}\t{group_name}\t{genes_num}\t{genes}\n')




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='汇总')
    parser.add_argument('--gene_table', help='基因统计表格')
    parser.add_argument('--taxonomy_tabel', help='物种注释表格')
    parser.add_argument('--result_suffix', help='输出结果前缀')

    args = vars(parser.parse_args())
    main(args)



