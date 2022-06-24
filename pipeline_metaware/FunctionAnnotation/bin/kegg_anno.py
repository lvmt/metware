#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-31 09:10:18
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-31 09:10:18
'''

'''
给KEGG比对结果添加注释信息 

'''


import pandas as pd  
import json 

def get_df(infile):
    return pd.read_csv(infile, sep='\t')



def main(args):
    m8_df = pd.read_csv(args['m8'], sep='\t', header=None).rename(columns = {
                                                                    0: 'Query', 
                                                                    1: 'Gene', 
                                                                    2: 'pident',
                                                                    3: 'length',
                                                                    4: 'mismatch',
                                                                    5: 'gapopen',
                                                                    6: 'qstart',
                                                                    7: 'qend',
                                                                    8: 'sstart',
                                                                    9: 'send',
                                                                    10: 'evalue',
                                                                    11: 'bitscore'}) 
    m8_df = m8_df.loc[:, ['Query', 'Gene']]  # 只需要指定行即可
    ko_gene_df = get_df(args['ko_gene'])
    m8_gene_com = pd.merge(m8_df, ko_gene_df, on='Gene', how='left').fillna('-')  # 合并m8和gene

    kegg_relation_info = json.load(open(args['ko_relation'])) 
    with open(args['result'], 'w') as fw:
        fw.write('Query\tGene\tKO\tKO_Name\tKO_Definition\tKO_EC\tModule\tKO_Pathway\n')
        for row in m8_gene_com.values:
            query, gene, ko = row 
            if ko == '-':
                ko_name = '-'
                ko_definition = '-'
                ko_ec = '-'
                ko_modual = '-'
                ko_pathway = '-'
            else:
                ko_name = kegg_relation_info[ko]['ko_name']
                ko_definition = kegg_relation_info[ko]['ko_definition']
                ko_ec = kegg_relation_info[ko]['ec']
                ko_modual = kegg_relation_info[ko]['module']
                ko_pathway = kegg_relation_info[ko]['pathway']

            fw.write(f'{query}\t{gene}\t{ko}\t{ko_name}\t{ko_definition}\t{ko_ec}\t{ko_modual}\t{ko_pathway}\n')




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='kegg 注释')
    parser.add_argument('--m8', help='序列比对结果文件')
    parser.add_argument('--ko_gene', help='kegg中基因与ko的对应关系-库文件')
    parser.add_argument('--ko_relation', help='kegg各个数据库间关系库文件')
    parser.add_argument('--result', help='注释输出结果文件')


    args = vars(parser.parse_args())
    main(args)

