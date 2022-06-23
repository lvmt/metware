#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-25 13:50:08
code check: 2022-06-22 成功推出诺禾计算逻辑
'''


# 注释结果及基因丰度表出发, 整合结果
import pandas as pd  
from collections import defaultdict


def get_group_info(sample_file):
    group_info = defaultdict(list)
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            sample = linelist[0]
            group = linelist[1]
            group_info[group].append(sample)
    return group_info


def get_group_df(df, group_info):
    for group in group_info:
        df[group] = df.loc[:, group_info[group]].sum(axis=1)
    
    df = df.loc[:, list(group_info.keys())]
    return df 


def main(args):
    group_info = get_group_info(args['sample_file'])
    tax_df = pd.read_csv(args['taxonomy_table'], sep='\t', header=None).rename(columns={0: 'Reference_ID', 1: 'taxonomy'})
    abun_df = pd.read_csv(args['abundance_table'], sep='\t')

    merge_df = pd.merge(tax_df, abun_df, on='Reference_ID', how='outer').fillna('k__Others;p__Others;c__Others;o__Others;f__Others;g__Others;s__Others')
    # 诺禾计算过程中的核心过程
    merge_df['taxonomy'] = merge_df['taxonomy'].apply(lambda x: x if len(x.split(';')) == 7 else ';'.join(x.split(';') + ['Other'] * (7 - len(x.split(';')))))
    # 某些基因在diamond比对过程中,会被舍弃掉，因此可能出现物种注释的基因和之前定量的基因数目存在差异
    # 在合并过程中,对于这些有差异的基因,认为其物种注释结果为others
    # 绝对不能将差异基因舍弃,否则前面的均一化结果就会失效,亲测
    all_result = '{result_suffix}.total.xls'.format(**args)
    merge_df.to_csv(all_result, sep='\t', index=None)
    class_info = {
        'kingdom': 1,
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'genus': 6,
        'species': 7
    }
    #################### 样本维度 进行层级拆分
    for _class in class_info:
        merge_df['Taxonomy'] = merge_df['taxonomy'].apply(lambda x:';'.join(x.split(';')[:class_info[_class]]))

        # 复现诺禾结果中, 对于Others的判断为, 该层级的为unclass或者others的归为Others
        # 核心代码
        merge_df.loc[merge_df['Taxonomy'].apply(lambda x: 'Unclassified' in x.split(';')[-1]), 'Taxonomy'] = 'Others'
        merge_df.loc[merge_df['Taxonomy'].apply(lambda x: 'Other' in x.split(';')[-1]), 'Taxonomy'] = 'Others'
        stat_df = merge_df.groupby('Taxonomy').agg('sum')
        
        ## 绝对丰都
        absolute_out = '{result_dir}/Absolute/{result_suffix}.absolute.{_class}.xls'.format(**args, **locals())
        absolute_df = stat_df.reset_index()
        absolute_df[_class] = absolute_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[-2:])) 
        absolute_df.to_csv(absolute_out, sep='\t', index=None) 
        
        ## 相对丰度
        relative_out = '{result_dir}/Relative/{result_suffix}.relative.{_class}.xls'.format(**args, **locals())
        relative_df = (stat_df / stat_df.sum(axis=0)).reset_index()
        relative_df[_class] = relative_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[-2:])) 
        relative_df.to_csv(relative_out, sep='\t', index=None)
        
        ## 分组维度 层级拆分（只进行相对丰度层面的）
        group_out = '{result_dir}/Relative/{result_suffix}.group.relative.{_class}.xls'.format(**args, **locals())
        group_df = stat_df / stat_df.sum(axis=0)
        group_df = get_group_df(group_df, group_info)
        group_df = (group_df / group_df.sum(axis=0)).reset_index()  # 重新归一化
        group_df[_class] = group_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[-2:])) 
        group_df.to_csv(group_out, sep='\t', index=None)




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='combine abundance and taxonomy')
    parser.add_argument('--abundance_table', help='基因丰度表格 Unigenes.readsNum.even.xls')
    parser.add_argument('--taxonomy_table', help='物种注释表格 Unigenes.lca.tax.xls')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--result_suffix', help='二者合并结果前缀')
    parser.add_argument('--result_dir', help='输出结果目录,会同时输出绝对和相对丰度结果')

    args = vars(parser.parse_args())
    main(args)