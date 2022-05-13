#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-11 11:28:40
'''


# MEGAN物种注释文件处理
import sys 

infile  = sys.argv[1]
outfile = sys.argv[2]


with open(infile, 'r') as fr, open(outfile, 'w') as fw:
    keys = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    for line in fr:
        tmp_item = {'k': 'other', 
                    'p': 'other', 
                    'c': 'other',
                    'o': 'other',
                    'f': 'other',
                    'g': 'other',
                    's': 'other'}

        linelist = line.strip('\n').split('\t')
        gene = linelist[0]
        taxonomy = linelist[1]
        tax_list = taxonomy.strip(';').split('; ')
        try:
            tax_list.remove('NCBI')
            tax_list.remove('cellular organisms')
        except Exception as e:
            pass 

        # 将unclass 转变为other
        for index in range(len(tax_list)):
            if 'unclass' in tax_list[index]:
                tax_list[index] = 'other'

        tmp_item.update(dict(zip(keys, tax_list)))

        new_taxonomy = ['{0}__{1}'.format(k,v) for k,v in tmp_item.items()]
        new_taxonomy = ';'.join(new_taxonomy)
        fw.write(f'{gene}\t{new_taxonomy}\n')


