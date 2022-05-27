#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-26 08:52:07
'''

# 基于丰度的绝对定量数据,绘制动态krona图

import pandas as pd  
import os  



def main(args):

    df = pd.read_csv(args['absolute_abundance_table'], sep='\t')
    df = df.drop(columns=['species'])
    samples = df.columns[1:]

    krona_file = []

    for sample in samples:
        tmpfile = '{out_dir}/{sample}.krona'.format(**args, **locals())
        krona_file.append(tmpfile)

        with open(tmpfile, 'w') as fw:
            sub_df = df.loc[:, ['Taxonomy', sample]]
            for row in sub_df.values:
                tax, gene_count = row 
                tax = tax.replace('k__', '').replace('p__','').replace('c__','').replace('o__','').replace('f__','').replace('g__','').replace('s__', '')
                # while tax.startswith('other') or tax.endswith('other'):
                #     tax = tax.strip(';other')
                # if not tax:
                #     continue

                tax = tax.replace('|', '\t')
                fw.write(f'{gene_count}\t{tax}\n')

    krona_files = ' '.join(krona_file)
    cmd = '/share/software/apps/anaconda3/envs/microbiome/bin/ktImportText {krona_files} -o {out_dir}/krona.html'.format(**locals(), **args)
    os.system(cmd)




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='绘制krona图')
    parser.add_argument('--absolute_abundance_table', help='输入, 使用种水平的文件')
    parser.add_argument('--out_dir', help='输出目录')

    args = vars(parser.parse_args())
    main(args)




