#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-20 09:57:56
'''


# 对预测得到的基因进行系列分析
'''
0. core-pan图
1. 基因数目箱线图
2. 韦恩图
3. 相关性图
'''


from matplotlib.transforms import Bbox
import pandas as pd  
from collections import defaultdict
import matplotlib.pyplot as plt 
plt.switch_backend('agg')



def get_df(infile):
    return pd.read_csv(infile, sep='\t', index_col=0)


def get_group_info(infile):
    group_info = defaultdict(list)
    with open(infile, 'r') as fr:
        for line in fr:
            linelist = line.strip('\n').split('\t')
            group = linelist[1]
            sample = linelist[0]
            group_info[group].append(sample)
    return group_info


def boxplot(df, group_info):
    # 分组箱线图
    '''
    group = {
        'A': [a1, a2, a3],
        'B': [b1, b2, b3]
    }
    '''
    gene_sum = defaultdict(list)
    print(group_info)
    for group in group_info:
        for sample in group_info[group]:
            gene_nums = sum(df[sample] > 0)  
            gene_sum[group].append(gene_nums)

    plot_df = pd.DataFrame(gene_sum)
    plot_df.boxplot().figure.savefig('hh.png', Bbox='tight', dpi=400)


def main(args):
    df = get_df(args['gene_table'])
    group_info = get_group_info(args['sample_file'])

    boxplot(df, group_info)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='绘图')
    parser.add_argument('--gene_table', help='基因汇总文件')
    parser.add_argument('--sample_file', help='配置文件')
    
    args = vars(parser.parse_args())

    main(args)


