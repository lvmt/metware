#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-26 14:34:42
@Last Modified by:   lvmengting
@Last Modified time: 2022-04-26 14:34:42
'''


'''
toolkits
'''

import os  
from collections import defaultdict



def mkdirs(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def write_cmd(cmd, shellname):
    dirname = os.path.dirname(shellname)
    mkdirs(dirname)
    with open(shellname, 'w') as fw:
        fw.write(cmd)



def get_fq_info(sample_file):
    '''
    sample  type fq1    fq2
    {
        'sample1': [fq1, fq2]
    }
    '''
    fq_info_dict = defaultdict(list)
    with open(sample_file, 'r') as fr:
        for line in fr:
            if not line.strip(''):
                continue
            linelist = line.strip('\n').split('\t')
            if linelist[0].lower() in ('sample'):
                continue
            fq_info_dict[linelist[0]] = [linelist[2], linelist[3]]
    
    return fq_info_dict


def get_group_info(sample_file):
    # 获取分析样本的属组

    pass 



    

