#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-12 16:34:11
'''

# 物种丰度前10绘制 
from matplotlib.font_manager import FontProperties
import pandas as pd  
import matplotlib.pyplot as plt 
plt.switch_backend('agg')



def top10(args):

    abundance_file = args.abundance_table
    figure_name = args.figure_name
    top = args.top

    # 数据库处理操作
    df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    top10_index = df.sum(axis=1).sort_values(ascending=False)[:top].index
    top10_df = df.loc[top10_index, :]
    other_df = pd.DataFrame(df.drop(labels=top10_index, axis=0).sum(axis=0)).T.rename(index={0:f'Notop{top}_Other'})
    plot_df = top10_df.append(other_df)
    plot_df_ratio = plot_df / plot_df.sum(axis=0)  # 转化为占比

    # # 处理索引
    new_index = {item:item.split(';')[-1] for item in plot_df_ratio.index}  # 简化图例
    plot_df_ratio = plot_df_ratio.rename(index=new_index)
    # print(plot_df_ratio)

    # 绘图
    plt.figure()
    p1 = plt.figure(figsize=(8, 6), dpi=80)
    plot_df_ratio = plot_df_ratio.T
    fig = plot_df_ratio.plot.bar(stacked=True, legend=True)
    plt.title('Taxonomy Abundance')
    plt.xlabel('Sample Name')
    plt.ylabel('Relative Abundance')
    plt.xticks(rotation=45)
    plt.legend(bbox_to_anchor=(0.80, -0.18), fontsize=8,title='taxonomy', ncol=2, frameon=False)
    fig.figure.savefig(figure_name, dpi=600, bbox_inches='tight')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='物种丰度top10绘图')
    parser.add_argument('--abundance_table', help='样本物种丰度统计表格')
    parser.add_argument('--figure_name', help='输出文件')
    parser.add_argument('--top', help='输出排名前多少的物种:默认top10', default=10, type=int)

    args = parser.parse_args()

    top10(args)




