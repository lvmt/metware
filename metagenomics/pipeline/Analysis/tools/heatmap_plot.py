#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-13 14:24:41
'''

# 物种基因数目热图 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import seaborn as sns
plt.switch_backend('agg')
# import palettable



def heatmap(args):
    # 丰度前35的物种
    # 绘制物种-基因数目热图
    # 绘制物种-丰度热图
    # sns.color_palette('hls', 8)  # 选取颜色数量
    abundance_file = args.abundance_table
    gene_file = args.gene_count_table 

    abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    gene_df = pd.read_csv(gene_file, sep='\t', index_col=0)
    plot_index = abundance_df.sum(axis=1).sort_values(ascending=False)[:10].index  # top物种的index
    abundance_df_plot = abundance_df.loc[plot_index, :]
    gene_df_plot = gene_df.loc[plot_index, :] 

    # 按照属水平进行颜色分类
    abundance_df_plot['class'] = [item.split(';')[-2] for item in abundance_df_plot.index] 
    gene_df_plot['class'] = [item.split(';')[-2] for item in gene_df_plot.index] 
    row_color = dict(zip(abundance_df_plot['class'].unique(), sns.color_palette('hls', len(abundance_df_plot['class'].unique()))))



    abundance_df_new_index = {item:item.split(';')[-1] for item in abundance_df_plot.index}
    abundance_df_plot = abundance_df_plot.rename(index=abundance_df_new_index)

    gene_df_plot_new_index = {item:item.split(';')[-1] for item in gene_df_plot.index}
    gene_df_plot = gene_df_plot.rename(index=gene_df_plot_new_index)

    # print(abundance_df_plot)

    plt.rcParams['font.sans-serif']=['SimHei']  # 用于显示中文
    plt.rcParams['axes.unicode_minus'] = False  # 用于显示中文

    plt.figure(dpi=200)
    abundance_fig = sns.clustermap(data=abundance_df_plot.iloc[:, :-1], col_colors=['r','r','r', 'r','g', 'g', 'g', 'g'], row_colors=abundance_df_plot['class'].map(row_color))
    abundance_fig.savefig('abundance_heatmap.png')

    gene_fig = sns.clustermap(data=gene_df_plot.iloc[:, :-1], col_colors=['r','r','r', 'r','g', 'g', 'g', 'g'], row_colors=abundance_df_plot['class'].map(row_color))
    gene_fig.savefig('gene_heatmap.png')







if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='热图绘制')
    parser.add_argument('--abundance_table', help='(属)水平物种丰度统计表')
    parser.add_argument('--gene_count_table', help='(属)水平物种基因统计表')
    parser.add_argument('--group_config', help='分组信息文件')
    parser.add_argument('--top', help='对排名前多少进行排序top30', default=30, type=int)
    parser.add_argument('--out-prefix', help='输出文件前缀')

    args = parser.parse_args()
    heatmap(args)


