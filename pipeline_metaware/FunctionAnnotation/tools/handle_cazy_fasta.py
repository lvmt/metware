#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-13 12:26:38
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-13 12:26:38
'''

# 基于protein.accession下载的cazy.fasta文件，部分fasta的名称无法对应上protein.accession里面的结果 
# 对fasta文件进行处理,只提取fasta名称和protein.accession一致的结果.


import Bio
from Bio import SeqIO
import argparse
import re  



parser = argparse.ArgumentParser(description='获取指定序列信息')
parser.add_argument('--access', help='accession numbers')
parser.add_argument('--fasta', help='fasta文件')
parser.add_argument('--result', help='结果文件')

args = vars(parser.parse_args())


def get_all_access(infile):
    access_names = []
    with open(infile, 'r') as fr:
        for line in fr:
            name = line.strip('\n')
            access_names.append(name)

    return access_names


def main(args):
    access_names = get_all_access(args['access'])

    with open(args['result'] + 'pass', 'w') as fp, \
        open(args['result'] + 'fail', 'w') as ff:
        
        for record in SeqIO.parse(args['fasta'], 'fasta'):
            name = record.name
            seq = record.seq 

            name = re.search(r'([A-Z0-9]+\.*\d*)', name).group()
            
            if name in access_names:
                fp.write(f'>{name}\n')
                fp.write(f'{seq}\n')
            else: 
                ff.write(f'{name}\n')
                ff.write(f'{seq}\n')



main(args)
