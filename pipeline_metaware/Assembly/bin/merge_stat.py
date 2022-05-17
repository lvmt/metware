#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-17 09:11:41
'''


# 合并全部样本的组装评估的结果文件

from cmath import inf
import pandas as pd  
from functools import reduce

from yaml import parse


def get_df(infile):
    return pd.read_csv(infile, sep='\t')


def main(args):
    dfs = []
    for infile in args['stat_files']:
        df = get_df(infile)
        dfs.append(df)

    combine_df = reduce(lambda x,y: pd.merge(x, y, on='Assembly', how='outer'), dfs)
    combine_df.to_csv(args['merge_result'], sep='\t', index=None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='合并组装评估结果')
    parser.add_argument('--stat_files', help='quast统计结果', nargs='+')
    parser.add_argument('--merge_result', help='输出结果')

    args = vars(parser.parse_args())
    main(args)

