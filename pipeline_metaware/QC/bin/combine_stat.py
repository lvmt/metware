#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-16 15:26:47
'''

# 合并全部样本的统计结果文件 

import argparse
import os
import pandas as pd

def get_df(infile):
    df = pd.read_csv(infile, sep='\t')
    name = os.path.basename(infile.replace('.stat.txt', ''))
    df['Sample'] = name
    df = df.loc[:,[df.columns.tolist()[-1]] + df.columns.tolist()[:-1]]
    return df


def main(args):
    df_list = []
    for infile in args['stat_files']:
        df = get_df(infile)
        df_list.append(df)

    merge_df = pd.concat(df_list, axis=0)
    merge_df.to_csv(args['stat_merge'], index=None, sep='\t')





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='合并统计文件')
    parser.add_argument('--stat_files', help='分析中使用到的配置文件', nargs='+')
    parser.add_argument('--stat_merge', help='全部样本的合并统计结果')

    args = vars(parser.parse_args())
    main(args)

