#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-25 13:50:08
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-25 13:50:08
'''


# 注释结果及基因丰度表出发, 整合结果
import pandas as pd  



def main(args):
    tax_df = pd.read_csv(args['taxonomy_table'], sep='\t')
    tax_df = tax_df.drop(columns=['taxid'])
    abun_df = pd.read_csv(args['abundance_table'], sep='\t')
    merge_df = pd.merge(tax_df, abun_df, on='gene', how='outer').fillna('k__Other|p__Other|c__Other|o__Other|f__Other|g__Other|s__Other')
    # 某些基因在diamond比对过程中,会被舍弃掉，因此可能出现物种注释的基因和之前定量的基因数目存在差异
    # 在合并过程中,对于这些有差异的基因,认为其物种注释结果为others
    # 绝对不能将差异基因舍弃,否则前面的均一化结果就会失效,亲测
    
    # 汇总结果 
    all_result = '{result_suffix}.total.xls'.format(**args)
    merge_df.to_csv(all_result, sep='\t', index=None)

    ## 进行层级拆分
    class_info = {
        'kingdom': 1,
        'phylum': 2,
        'class': 3,
        'order': 4,
        'family': 5,
        'Genus': 6,
        'species': 7
    }
    for _class in class_info:
        merge_df['Taxonomy'] = merge_df['taxonomy'].apply(lambda x:'|'.join(x.split('|')[:class_info[_class]]))
        class_df = merge_df.groupby('Taxonomy').agg('sum').reset_index()
        class_df[_class] = class_df['Taxonomy'].apply(lambda x:'|'.join(x.split('|')[-2:])) 

        out = '{result_suffix}.{_class}.xls'.format(**args, **locals())
        class_df.to_csv(out, sep='\t', index=None)
        



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='combine abundance and taxonomy')
    parser.add_argument('--abundance_table', help='基因丰度表格')
    parser.add_argument('--taxonomy_table', help='物种注释表格')
    parser.add_argument('--result_suffix', help='二者合并结果前缀')

    args = vars(parser.parse_args())
    main(args)