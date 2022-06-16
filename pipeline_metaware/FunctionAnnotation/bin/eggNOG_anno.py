#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-15 14:57:01
'''

# eggNOG 数据库注释 
import json  


def main(args):
    relation_info = json.load(open(args['eggNOG_relation']))
    with open(args['m8'], 'r') as fr, open(args['result'], 'w') as fw:
        fw.write('query\tgene\tOrtholog_Group\tFunctional_Category\tOG_Description\n')
        for line in fr:
            linelist = line.strip('\n').split('\t')
            query = linelist[0]
            gene = linelist[1]
            ogs = relation_info.get(gene)
            if ogs:
                for og in ogs:
                    OG, taxid, OG_class, OG_desc = og.split('|')
                    fw.write(f'{query}\t{gene}\t{OG}\t{OG_class}\t{OG_desc}\n')

      


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='eggNOG数据库注释')
    parser.add_argument('--m8', help='diamond cazy 比对结果')
    parser.add_argument('--eggNOG_relation', help='eggNOG关系数据库')
    parser.add_argument('--result', help='eggNOG注释结果')

    args = vars(parser.parse_args())

    main(args)