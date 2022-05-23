#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-23 10:47:13
'''


# core-pan基因稀释曲线
'''
Core gene是所有样本中均存在的同源基因
Pan gene泛基因包含所有非共有与共有基因
绘制结果的横坐标为样本数目, 因此在样本数目的维度上,需要进行组合
当样本总数过多时,生成的组合类型将会非常多;因此在绘制图片时,限制最大组合数为20.
'''


import matplotlib.pyplot  as plt
from matplotlib.transforms import Bbox
import pandas as pd  
import random
from scipy.special import comb
from collections import defaultdict
import seaborn as sns 
plt.switch_backend('agg')


def sample_combine_info(sample_file):
    '''
    样本组合关系
    {
        1: ['HC1', 'HK1'],
        2: ['HC1|HK1', 'HH1|HC2'],
    }
    '''
    sample_info = {}
    samples = []
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            samples.append(linelist[0])

    for i in range(1, len(samples) + 1):
        sample_info[i] = set()
        total_nums = comb(len(samples), i)
        need_nums = min(total_nums, 20)
        while True:
            tmp = '|'.join(sorted(random.sample(samples, i)))
            sample_info[i].add(tmp)
            if len(sample_info[i]) == need_nums:
                break

    return sample_info
    

def handle(gene_table, sample_info):
    '''
    nums_stat = {
        1: [nums, nums2],
        2: [nums, nums]
    }

    geneid_stat = {
        1(HC1): 'gene1|gene2',
        1(HC2): 'gene1|gene2',
        2(HC1|HK1): 'gene1|gene2' 
    }
    '''
    core_nums_stat = defaultdict(list)
    pan_nums_stat = defaultdict(list)
    core_gene_file = '{result_suffix}.core_geneid'.format(**args)
    pan_gene_file = '{result_suffix}.pan_geneid'.format(**args)
    core_nums_plot = '{result_suffix}.core.gene.png'.format(**args)
    pan_nums_plot = '{result_suffix}.pan.gene.png'.format(**args)

    df = pd.read_csv(gene_table, sep='\t')
    with open(core_gene_file, 'w') as fw_core, open(pan_gene_file, 'w') as fw_pan:
        for sample_num in sample_info:  # sample_num: 1/2/3...
            for item in sample_info[sample_num]:
                samples = item.split('|')  # item: 'HC1|HK1'
                sub_df = df.loc[:, samples]
                sub_df = sub_df > 0
                core_index = sub_df[sub_df.sum(axis=1) == len(sub_df.columns)].index  # 共有基因才输出
                pan_index = sub_df[sub_df.sum(axis=1) > 0].index    # 有就输出

                core_nums_stat[sample_num].append(len(core_index))
                pan_nums_stat[sample_num].append(len(pan_index))
                core_geneid = '|'.join(df.loc[core_index, ]['gene'].values)
                pan_geneid = '|'.join(df.loc[pan_index, ]['gene'].values)
                fw_core.write(f'{item}\t{core_geneid}\n')
                fw_pan.write(f'{item}\t{pan_geneid}\n')

    # 绘图
    core_df = pd.DataFrame.from_dict(core_nums_stat, orient='index').T
    pan_df = pd.DataFrame.from_dict(pan_nums_stat, orient='index').T
    sns.boxplot(data=core_df)
    plt.xlabel('Number of samples')
    plt.ylabel('Number of non-redundant genes')
    plt.savefig(core_nums_plot, dpi=400, Bbox='tight')
    sns.boxplot(data=pan_df)
    plt.xlabel('Number of samples')
    plt.ylabel('Number of non-redundant genes')
    plt.savefig(pan_nums_plot, dpi=400, Bbox='tight')





def main(args):
    sample_info = sample_combine_info(args['sample_file'])
    handle(args['gene_table'], sample_info)


     





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='core-pan基因统计分析')
    parser.add_argument('--sample_file', help='输入: 样本信息')
    parser.add_argument('--gene_table', help='输入: 非冗余基因列表')
    parser.add_argument('--result_suffix', help='输出文件前缀')

    args = vars(parser.parse_args())
    main(args)



