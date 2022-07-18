#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-07-14 08:34:02
'''
from collections import defaultdict


def file2dict(infile):
    info = defaultdict(list)
    with open(infile, 'r') as fr:
        for line in fr:
            linelist = line.strip('\n').split('\t')
            query = linelist[0]
            info[query].append(linelist)
    return info 


def get_filter_result(_list):
    # sourcery skip: inline-immediately-returned-variable, list-comprehension, move-assign-in-block
    # _list 为一个query的全部比对结果
    filter_list = []
    evalue_list = [item[-2] for item in _list]
    min_evalue = min(map(float, evalue_list)) * 10 # 
    
    for item in _list:
        if float(item[-2]) <= min_evalue:
            filter_list.append(item)
    return filter_list
    



def main(args):
    info = file2dict(args['m8'])
    
    with open(args['filter_result'], 'w') as fw:
        for query in info:
            _list = info[query]
            filter_list = get_filter_result(_list)
            
            for row in filter_list:
                fw.write('{}\n'.format('\t'.join(row)))



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='比对结果基于最小evalue值进行过滤')
    parser.add_argument('--m8', help='输入文件')
    parser.add_argument('--filter_result', help='输出文件, 对m8文件进行过滤后的结果')

    args = vars(parser.parse_args())
    main(args)


        
