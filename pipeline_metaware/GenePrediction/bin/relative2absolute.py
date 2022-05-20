#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-19 14:25:22
'''


# 将相对基因的相对丰度结果转换为绝对丰度结果
'''
转换思想:
统计基因在各样品中的覆盖 reads 数表中每列基因数之和，选取最大值作为放大倍数，
以此最大值×基因在各样品中的相对丰度表 Unigenes.readsNum.relative.xls 得出的相对丰度即
可得到均一化后基因的绝对丰度。
大概意思就是: 将相对丰度值整体乘以一个系数
'''


import pandas as pd  



def main(args):
    max_coeff = pd.read_csv(args['gene_table'], sep='\t', index_col=0).sum(axis=0).max()
    df = pd.read_csv(args['relative_table'], sep='\t', index_col=0)
    df = df * max_coeff
    df = df.reset_index()
    df.to_csv(args['absolute_table'], sep='\t', index = None)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='相对丰度转化为绝对丰度')
    parser.add_argument('--relative_table', help='输入文件,相对定量结果')
    parser.add_argument('--gene_table', help='输入文件, 基因计数结果')
    parser.add_argument('--absolute_table', help='输出结果, 绝对定量结果')

    args = vars(parser.parse_args())
    main(args)






