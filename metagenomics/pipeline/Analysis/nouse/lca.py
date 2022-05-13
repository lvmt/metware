#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-10 13:54:21
'''


# 借助taxonkit进行LCA分析
from ast import arg
from ctypes import alignment
from operator import ge
import pandas as pd
from collections import defaultdict
import subprocess

from zmq import DEALER


def get_alignment_df(align):
    df = pd.read_csv(align, sep='\t', header=None, dtype=str, names=list(range(12)))
    df = df.loc[:, [0,1]].rename(columns={0:'gene', 1: 'access'})
    return df  


def get_relation_df(relation):
    return pd.read_csv(relation, dtype=str, sep='\t').rename(columns={'accession.version': 'access'})


def get_merge_df(align_df, relation_df):
    return pd.merge(align_df, relation_df, on='access', how='outer')


def get_gene_taxid(merge_df):
    gene_taxid = defaultdict(list)

    for row in merge_df.values:
        gene, acc, taxid = row 
        gene_taxid[gene].append(taxid)

    return gene_taxid


def lca_analysis(gene_taxid):
    gene_taxid_only = {}
    for gene in gene_taxid:
        taxids = ','.join(gene_taxid[gene])
        cmd = f'echo {taxids}| ~/.conda/envs/python3_lmt/bin/taxonkit lca -s ","'
        final_taxid = subprocess.getoutput(cmd).split('\t')[-1]
        gene_taxid_only[gene] = final_taxid

    return gene_taxid_only


def main(args):
    alignment_df = get_alignment_df(args.alignment)
    relation_df = get_relation_df(args.relation)
    merge_df = get_merge_df(alignment_df, relation_df)
    gene_taxid = get_gene_taxid(merge_df)
    final_gene_taxid = lca_analysis(gene_taxid)

    with open(args.result, 'w') as fw:
        fw.write('gene\ttaxid\n')
        for gene in final_gene_taxid:
            taxid = final_gene_taxid[gene]
            fw.write(f'{gene}\t{taxid}\n')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='LCA analysis')
    parser.add_argument('--alignment', help='diamond比对结果,m8')
    parser.add_argument('--relation', help='accession与taxid关系文件')
    parser.add_argument('--result', help='结果文件,一条序列对应唯一taxid')

    args = parser.parse_args()

    main(args)

    






