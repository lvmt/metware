#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-06 08:35:57
'''


from itertools import count
from matplotlib.pyplot import legend
import pandas as pd
from requests import head 



def get_count_df(count_file):
    return pd.read_csv(count_file, sep=',', header=None, names=['gene', 'count'])


def get_length_df(length_file):
    return pd.read_csv(length_file, sep=',', header=None, names=['gene', 'length'])


def main(count_df, length_df, outfile):
    merge_df = pd.merge(count_df, length_df, on='gene', how='inner')  # 丢掉基因count不满足阈值的基因
    merge_df['count'] = merge_df['count'].values.astype(float)
    merge_df['length'] = merge_df['length'].values.astype(float)

    # 归一化分母
    merge_df['average'] = sum(merge_df.apply(lambda row: row['count'] / row['length'], axis=1))
    merge_df['abundance'] = merge_df.apply(lambda row: float(row['count']) / row['length'] / row['average'], axis=1)

    merge_df.to_csv(outfile, sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='丰度统计')
    parser.add_argument('--gene_count', help='基因reads计数文件')
    parser.add_argument('--gene_length', help='基因长度文件')
    parser.add_argument('--abundance_table', help='输出丰度统计结果')

    args = vars(parser.parse_args())
    
    count_df = get_count_df(args['gene_count'])
    length_df = get_length_df(args['gene_length'])
    main(count_df, length_df, args['abundance_table'])









