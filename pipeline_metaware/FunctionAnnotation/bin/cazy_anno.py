#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-13 11:28:14
'''

# 碳水化合物数据库注释 
import json  




def main(args):
    relation_info = json.load(open(args['cazy_relation']))
    with open(args['m8'], 'r') as fr, open(args['result'], 'w') as fw:
        fw.write('query\tgene\tfamily\tec\n')
        for line in fr:
            linelist = line.strip('\n').split('\t')
            query = linelist[0]
            gene = linelist[1]
            family = relation_info[gene]['family']
            ec = relation_info[gene]['ec']

            fw.write('{query}\t{gene}\t{family}\t{ec}\n'.format(**locals()))



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='碳水化合物注释')
    parser.add_argument('--m8', help='diamond cazy 比对结果')
    parser.add_argument('--cazy_relation', help='cazy关系数据库')
    parser.add_argument('--result', help='cazy注释结果')

    args = vars(parser.parse_args())

    main(args)
    


