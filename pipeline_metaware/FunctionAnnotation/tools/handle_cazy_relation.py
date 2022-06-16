#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-13 10:55:12
'''

'''
制作cazy关系数据库,用于anno
{
    'accession':
        'family': [f1, f2], # 可能对应多个家族
        'ec': [ec1, ec2],

}
'''

from  collections import defaultdict
from hashlib import new
import json

import argparse
parser = argparse.ArgumentParser(description='处理cazy数据库文件')
parser.add_argument('--family', help='family relation')
parser.add_argument('--ec', help='ec关系文件')
parser.add_argument('--result', help='输出结果')

args = vars(parser.parse_args())


info = defaultdict(set)
with open(args['family'], 'r') as fr:
    for line in fr:
        linelit = line.strip('\n').split('\t')
        fam = linelit[0]
        acc = linelit[1]

        if acc in info:
            info[acc]['family'].add(fam)
        else:
            info[acc]['family'] = set()
            info[acc]['ec'] = set()

            info[acc]['family'].add(fam)


with open(args['ec'], 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        acc = linelist[0]
        ec = linelist[1]

        if acc in info:
            info[acc]['ec'].add(ec)
        else:
            info[acc]['family'] = set()
            info[acc]['ec'] = set()
            info[acc]['ec'].add(ec)

# 改写info
new_info = defaultdict(dict)
for acc in info:
    fam = '|'.join(info[acc]['family'])
    ec = '|'.join(info[acc]['ec'])

    new_info[acc]['family'] = fam
    new_info[acc]['ec'] = ec 


with open(args['result'], 'w', encoding='utf-8') as fw:
    json.dump(new_info, fw) 


    


