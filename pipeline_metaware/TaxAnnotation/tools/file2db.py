#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-07-14 09:54:44
'''


# 将物种注释文件存入数据库中


import sqlite3
import argparse



parser = argparse.ArgumentParser(description='物种文件存入输入库中')
parser.add_argument('--db', help='连接的数据库')
parser.add_argument('--infile', help='待写入数据库的文件')


args = vars(parser.parse_args())

con = sqlite3.connect(args['db'])
cursor = con.cursor()



with open(args['infile'], 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        acc = linelist[0]
        taxid = linelist[1]
        
        cursor.execute(f"INSERT INTO RELATION VALUES ('{acc}', '{taxid}');")
        con.commit()





