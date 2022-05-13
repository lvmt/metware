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

def mkdirs(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def write_cmd(cmd, shellname):
    dirname = os.path.dirname(shellname)
    mkdirs(dirname)
    with open(shellname, 'w') as fw:
        fw.write(cmd)


    
