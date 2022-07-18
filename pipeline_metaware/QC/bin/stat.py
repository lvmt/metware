#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-16 14:53:21
'''


# 统计原始数据与最终质控后的数据结果


import pandas as pd 


def get_info(infile):
    return pd.read_csv(infile, sep='\t').loc[0,].to_dict()


def main(args):

    raw1_info = get_info(args['raw1'])
    raw2_info = get_info(args['raw2'])
    clean1_info = get_info(args['clean1'])
    clean2_info = get_info(args['clean2'])


    raw_total_bases = raw1_info['bases'] + raw2_info['bases']
    row_total_Q20 = round((raw1_info['Q20'] + raw2_info['Q20']) / raw_total_bases * 100, 2)
    row_total_Q30 = round((raw1_info['Q30'] + raw2_info['Q30']) / raw_total_bases * 100, 2)
    row_total_GC = round((raw1_info['GC'] + raw2_info['GC']) / raw_total_bases * 100, 2)
    raw_total_bases = round(raw_total_bases / 1000 / 1000, 2)

    clean_total_bases = clean1_info['bases'] + clean2_info['bases']
    clean_total_Q20 = round((clean1_info['Q20'] + clean2_info['Q20']) / clean_total_bases * 100, 2)
    clean_total_Q30 = round((clean1_info['Q30'] + clean2_info['Q30']) / clean_total_bases * 100, 2)
    clean_total_GC = round((clean1_info['GC'] + clean2_info['GC']) / clean_total_bases * 100, 2)
    clean_total_bases = round(clean_total_bases / 1000 / 1000, 2)

    with open(args['stat_result'], 'w') as fw:
        fw.write('RawData_bases(M)\tRawData_Q20\tRawData_Q30\tRawData_GC\tCleanData_bases(M)\tCleanData_Q20\tCleanData_Q30\tCleanData_GC\n')
        fw.write(f'{raw_total_bases}\t{row_total_Q20}\t{row_total_Q30}\t{row_total_GC}\t{clean_total_bases}\t{clean_total_Q20}\t{clean_total_Q30}\t{clean_total_GC}\n')





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='质控信息统计')

    parser.add_argument('--raw1', help='rawdata read1统计信息')
    parser.add_argument('--raw2')
    parser.add_argument('--clean1')
    parser.add_argument('--clean2')
    parser.add_argument('--stat_result')

    args = vars(parser.parse_args())
    main(args)

