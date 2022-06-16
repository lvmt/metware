#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-30 14:56:30
'''


# 构建关系表: 主键为KO
'''
    1 KO
    2 KO_name
    3 KO_description
    4 pathway_level1
    5 pathway_level2 
    6 pathway_level3(pathway-kid)
    7 geneid
    8 module
    9 enmyze
'''

import re  
import json 
from collections import defaultdict
from pprint import pprint


# with open('KO', 'r') as fr, open('ko.simple', 'w') as fw:
#     for line in fr:
#         linelist = line.strip('\n').split('\t')
#         ko = linelist[0]
#         info_list = re.split(r':|;|\[', linelist[-1])
#         name = info_list[0]
#         try:
#             dd = info_list[1]
#         except Exception as e:
#             dd = info_list[-1]

#         fw.write('{ko}\t{name}\t{dd}\n'.format(**locals()))



# 将kegg 库文件存储为json格式 
info = defaultdict(dict)
with open('ko_pathway_level.ec.mo.ko_name.xls', 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        ko = linelist[0]
        ko_name = linelist[1]
        ko_desc = linelist[2]
        pathway = linelist[3]
        level1 = linelist[4]
        level2 = linelist[5]
        level3 = linelist[6]
        ec = linelist[7]
        modual = linelist[8]

        if ko in info:
            info[ko]['ko_name'].add(ko_name)
            info[ko]['ko_definition'].add(ko_desc)
            info[ko]['pathway'].add(f'{pathway};{level1};{level2};{level3}')
            info[ko]['ec'].add(ec)
            info[ko]['modual'].add(modual)

        else:
            info[ko]['ko_name'] = set()
            info[ko]['ko_definition'] = set()
            info[ko]['pathway'] = set()
            info[ko]['ec'] = set()
            info[ko]['modual'] = set()
            info[ko]['ko_name'].add(ko_name)
            info[ko]['ko_definition'].add(ko_desc)
            info[ko]['pathway'].add(f'{pathway};{level1};{level2};{level3}')
            info[ko]['ec'].add(ec)
            info[ko]['modual'].add(modual)

new_info = defaultdict(dict)

for ko in info:
    new_info[ko]['ko_name'] = '|'.join(info[ko]['ko_name'])
    new_info[ko]['ko_definition'] = '|'.join(info[ko]['ko_definition'])
    new_info[ko]['pathway'] = '|'.join(info[ko]['pathway'])
    new_info[ko]['ec'] = '|'.join(info[ko]['ec'])
    new_info[ko]['modual'] = '|'.join(info[ko]['modual'])

with open('kegg_relation.json', 'w', encoding='utf-8') as f:
    json.dump(new_info, f)
pprint(new_info)