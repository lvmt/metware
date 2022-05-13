#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-10 14:59:57
'''

import pandas as pd   


def main(args):
    taxonomy_df = pd.read_csv(args.taxid_taxonomy, dtype=str, sep='\t')
    gene_df = pd.read_csv(args.gene_taxid, dtype=str, sep='\t')  # gene taxid
    taxonomy_df['taxonomy'] = taxonomy_df['kindom'].str.cat((taxonomy_df['phylum'], taxonomy_df['class'], taxonomy_df['order'], taxonomy_df['family'], taxonomy_df['genus'], taxonomy_df['species']),sep=';')
    taxonomy_df = taxonomy_df.loc[:, ['taxid', 'taxonomy']]

    pd.merge(gene_df, taxonomy_df, on='taxid', how='outer').to_csv(args.taxonomy_anno_table, sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='a')
    parser.add_argument('--gene_taxid', help='gene与taxid对应关系表')
    parser.add_argument('--taxid_taxonomy', help='taxid与物种注释关系表格')
    parser.add_argument('--taxonomy_anno_table', help='结果文件,物种注释表格')

    args = parser.parse_args()

    main(args)
    