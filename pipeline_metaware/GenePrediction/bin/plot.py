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

import random 
from functools import reduce
from matplotlib.transforms import Bbox
import pandas as pd
from scipy import rand  
import seaborn as sns 
from collections import defaultdict
import matplotlib.pyplot as plt 
plt.switch_backend('agg')



def get_samples(sample_file):
    samples = []
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            samples.append(linelist[0])
    return samples 


def get_group_info(infile):
    group_info = defaultdict(list)
    with open(infile, 'r') as fr:
        for line in fr:
            linelist = line.strip('\n').split('\t')
            group = linelist[1]
            sample = linelist[0]
            group_info[group].append(sample)
    return group_info


def get_target_group_df(df, group_info):
    # 获取分组df
    for _group in group_info:
        df[_group] = df.loc[:, group_info[_group].sum(axis=1)]

    df = df.loc[:, list(group_info.keys())]
    return df 


def boxplot(df, group_info, result_suffix):
    # 分组箱线图
    '''
    group = {
        'A': [a1, a2, a3],
        'B': [b1, b2, b3]
    }
    '''
    gene_sum = defaultdict(list)
    for group in group_info:
        for sample in group_info[group]:
            gene_nums = sum(df[sample] > 0)  
            gene_sum[group].append(gene_nums)

    boxplot_result = f'{result_suffix}.boxplot.png'
    plt.figure(figsize=(10, 10), dpi=100)
    ax = sns.boxplot(data=pd.DataFrame(gene_sum))
    ax.set_xlabel('Group')
    ax.set_ylabel('Count')
    ax.get_figure().savefig(boxplot_result, Bbox='tight')


def core_pan(df, samples, result_suffix):
    '''
    # core 核心基因, 样本共有基因(基因交集)
    # pan 泛基因, 样本基因并集 
    stat_info = {
        1: [nums1, nums2, ...],
        2: []
    }
    '''
    core_stat_info = defaultdict(list)
    pan_stat_info = defaultdict(list)
    max_random_nums = max(len(df.columns), 100)  # 设置随机抽样次数
    # gene_name_list = [list(df[df[sample] > 0].index) for sample in samples]
    # 耗时步骤
    # for sample_nums in range(1, len(samples) + 1):  # 随机抽取多少个样本
    #     for _ in range(max_random_nums):
    #         tmp = random.sample(gene_name_list, sample_nums)
    #         core_stat_info[sample_nums].append(len(reduce(lambda x,y: set(x) & set(y), tmp)))
    #         pan_stat_info[sample_nums].append(len(reduce(lambda x,y: set(x) | set(y), tmp)))
    random.seed(10)  # 保证结果可重现性
    for sample_nums in range(1, len(samples) + 1):
        for _ in range(max_random_nums):
            core_num = sum(df[random.sample(samples, sample_nums)].sum(axis=1) == sample_nums)
            pan_num = sum(df[random.sample(samples, sample_nums)].sum(axis=1) > 0)
            core_stat_info[sample_nums].append(core_num)
            pan_stat_info[sample_nums].append(pan_num)
            
    fig_core, ax_core = plt.subplots(1, 1, figsize=(20, 10))
    sns.boxplot(data=pd.DataFrame.from_dict(core_stat_info), ax=ax_core)
    ax_core.set_xlabel('Sample Numbers')
    ax_core.set_ylabel('Count')
    ax_core.set_title('Core gene analysis')
    fig_core.savefig(f'{result_suffix}.core.png', Bbox='tight', dpi=400)

    fig_pan, ax_pan = plt.subplots(1, 1, figsize=(20, 10))
    sns.boxplot(data=pd.DataFrame.from_dict(pan_stat_info), ax=ax_pan)
    ax_pan.set_xlabel('Sample Numbers')
    ax_pan.set_ylabel('Count')
    ax_pan.set_title('Pan gene analysis')
    fig_pan.savefig(f'{result_suffix}.pan.png', Bbox='tight', dpi=400) 



def main(args):
    df = pd.read_csv(args['gene_table'], sep='\t', index_col=0)
    samples = get_samples(args['sample_file'])
    group_info = get_group_info(args['sample_file'])

    # 将基因数大于0的设为1，方便后续统计
    for sample in samples:
        df.loc[df.loc[:, sample] > 0, sample] = 1
    
    # 分组基因数箱线图
    boxplot(df, group_info, args['result_suffix'])
    core_pan(df, samples, args['result_suffix'])




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='绘图')
    parser.add_argument('--gene_table', help='基因汇总文件')
    parser.add_argument('--sample_file', help='配置文件')
    parser.add_argument('--result_suffix', help='输出结果前缀')
    
    args = vars(parser.parse_args())

    main(args)


