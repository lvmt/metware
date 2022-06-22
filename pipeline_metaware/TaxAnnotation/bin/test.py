#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-22 11:04:17
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-22 11:04:17
'''

# 测试诺禾在物种拆分时，是如何处理others
import pandas as pd 
others = ['Others'] 

# def main(args):
#     with open(args['infile'], 'r') as fr, open(args['out'], 'w') as fw:
#         for line in fr:
#             if line.startswith('Reference_ID'):
#                 fw.write(line)
#                 continue
#             linelist = line.strip('\n').split('\t')
#             tax_list = linelist[13].split(';')
#             new_tax_list = tax_list + others * (7 - len(tax_list))
#             new_tax_list = ';'.join(new_tax_list)
#             linelist[13] = new_tax_list
#             fw.write('{}\n'.format('\t'.join(linelist)))


def stat(infile):
    class_info = {
        'kingdom': 1,
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'Genus': 6,
        'species': 7
    }

    merge_df = pd.read_csv(infile, sep='\t')
    merge_df['Taxonomy'] = merge_df['Taxonomy'].apply(lambda x: x if len(x.split(';')) == 7 else ';'.join(x.split(';') + ['Others'] * (7 - len(x.split(';')))))

    for _class in class_info:
        out = f'{infile}.{_class}.xls'
        merge_df['taxonomy'] = merge_df['Taxonomy'].apply(lambda x:';'.join(x.split(';')[:class_info[_class]]))

        merge_df.loc[merge_df['taxonomy'].apply(lambda x: 'Unclassified' in x.split(';')[-1]), 'taxonomy'] = 'Others'
        merge_df.loc[merge_df['taxonomy'].apply(lambda x: 'Others' in x.split(';')[-1]), 'taxonomy'] = 'Others'
        stat_df = merge_df.groupby('taxonomy').sum().reset_index() 
        stat_df[_class] = stat_df['taxonomy'].apply(lambda x: ';'.join(x.split(';')[-2:]))

        stat_df.to_csv(out, sep='\t', index=None)


stat('Unigenes.absolute.total.tax.xls')






