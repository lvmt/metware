#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-25 09:54:40
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-25 09:54:40
'''


'''
从库文件中taxonomy中提取制定taxid的层级注释信息
'''

from hamcrest import none
import pandas as pd  



def main(args):
    db = pd.read_csv(args['taxonomy_database'], sep='\t', header=None).rename(columns={0: 'taxid', 1: 'taxonomy1', 2: 'taxonomy'})
    db = db.drop(columns=['taxonomy1'])
    lca = pd.read_csv(args['lca'], sep='\t', header=None).rename(columns={0: 'gene', 1: 'taxid'})

    df = pd.merge(db, lca, on='taxid', how='right').fillna('k__Other|p__Other|c__Other|o__Other|f__Other|g__Other|s__Other')
    df.to_csv(args['result'], sep='\t', index=None)




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='taxid层级注释')
    parser.add_argument('--lca', help='物种注释软件生成的结果')
    parser.add_argument('--taxonomy_database', help='人为构建好的taxid-taxonomy数据库')
    parser.add_argument('--result', help='输出结果')

    args = vars(parser.parse_args())

    main(args)
