#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-07-13 16:20:43
@Last Modified by:   lvmengting
@Last Modified time: 2022-07-13 16:20:43
'''

# LCA 核心：树中两个节点的最低公共节点

import sys  
from collections import defaultdict


infile = sys.argv[1]
outfile = sys.argv[2]



def lca(_list):
    lca_item = []
    if len(_list) == 1:
        return ';'.join(_list[0])

    for index in range(min(map(len, _list))):
        if len({item[index] for item in _list}) == 1:
            lca_item.append(_list[0][index])
        else:
            break 

    return ';'.join(lca_item) if lca_item else 'Unclassified'



info = defaultdict(list)
with open(infile, 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        query = linelist[0]
        tax = linelist[-1].split(';')
        info[query].append(tax)


with open(outfile, 'w') as fw:
    fw.write('gene\ttaxonomy\n')
    for query in info:
        tax = lca(info[query])
        fw.write(f'{query}\t{tax}\n')


