#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-12 14:35:11
'''


# 绘制krona动态图
import pandas as pd  
import os  


def plot_krona(args):
    infile = args.infile
    outdir = args.outdir

    df = pd.read_csv(infile, sep='\t').fillna('0')
    samples = df.columns[1:]

    krona_file = []

    for sample in samples:
        tmpfile = f'{outdir}/{sample}.krona'
        krona_file.append(tmpfile)

        with open(tmpfile, 'w') as fw:
            sub_df = df.loc[:, ['taxonomy', sample]]
            for row in sub_df.values:
                tax, gene_count = row 
                tax = tax.replace('k__', '').replace('p__','').replace('c__','').replace('o__','').replace('f__','').replace('g__','').replace('s__', '')
                while tax.startswith('other') or tax.endswith('other'):
                    tax = tax.strip(';other')
                if not tax:
                    continue

                tax = tax.replace(';', '\t')
                fw.write(f'{gene_count}\t{tax}\n')

    krona_files = ' '.join(krona_file)
    cmd = f'/share/software/apps/anaconda3/envs/microbiome/bin/ktImportText {krona_files} -o {outdir}/krona.html'
    os.system(cmd)
    


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='绘制krona动态图')
    parser.add_argument('--infile', help='输入文件:taxonomy.geneNums.table')
    parser.add_argument('--outdir', help='输出目录')

    args = parser.parse_args()

    plot_krona(args)
