#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-07-15 10:51:38
'''


# 基于access - taxid - taxonomy间的关系，给query基因增加物种注释信息 

import pandas as pd   
import argparse



def main(args):
    
    m8_df = pd.read_csv(args['m8_filter'], sep='\t', header=None).rename(columns={1: 'access'})
    acc2taxid_df = pd.read_csv(args['acc2taxid'], sep='\t', header=None).rename(columns={0: 'access', 1: 'taxid'})
    merge_m8_taxid = pd.merge(m8_df, acc2taxid_df, on='access')
    
    tax_df = pd.read_csv(args['taxid2taxonomy'], sep='\t', header=None).rename(columns={0: 'taxid', 1: 'taxonomy'})
    merge_final = pd.merge(merge_m8_taxid, tax_df, on='taxid')
    
    merge_final.sort_values(by=0).to_csv(args['m8_taxonomy'], sep='\t', index=None, header=None)
    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='物种注释')
    parser.add_argument('--m8_filter', help='经过过滤的alignment文件')
    parser.add_argument('--acc2taxid', help='上一步骤求取得到的关系文件')
    parser.add_argument('--taxid2taxonomy', help='物种id与物种注释关系的配置文件')
    parser.add_argument('--m8_taxonomy', help='输出结果, m8文件增加注释信息的结果文件')

    args = vars(parser.parse_args())

    main(args)

