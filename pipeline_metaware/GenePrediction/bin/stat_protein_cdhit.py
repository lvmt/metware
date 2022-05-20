#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-19 10:14:37
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-19 10:14:37
'''


# 对蛋白去冗余的结果进行统计分析
'''
rep_id  Len(aa) Gene_Nums   Gene_ids
'''

import re 


def main(args):
    infile = args['clstr']
    outfle =args['stat']
    info  = {}
    with open(infile, 'r') as fr:
        for line in fr:
            line = line.strip('\n')
            linelist = line.strip('\n').split('\t')
            if line.startswith('>'):
                name = linelist[0]
                info[name] = []
            else:
                tmp = linelist[-1]
                length = re.search(r'(\d+)aa,', tmp).group(1)
                gene_name = re.search(r'(>.+)\.\.\.', tmp).group(1)

                string = f'{length};{gene_name}'
                if line.endswith('*'):
                    info[name].insert(0, string)
                else:
                    info[name].append(string)      

    with open(outfle, 'w') as fw:
        fw.write('Rep_id\tLen(nt/aa)\tNum\tSeq_ID\n')
        for clst in info:
            gene_ids = info[clst]
            rep = gene_ids[0].split(';')[1]  # 代表性基因
            rep_len = gene_ids[0].split(';')[0]
            gene_nums = len(gene_ids)
            gene_ids = ','.join([gene.split(';')[1] for gene in gene_ids])
            
            fw.write(f'{rep}\t{rep_len}\t{gene_nums}\t{gene_ids}\n')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='CD_HIT(protein level) 统计分析')
    parser.add_argument('--clstr', help='去冗余过程中生成的聚类文件')
    parser.add_argument('--stat', help='输出的统计结果文件')

    args = vars(parser.parse_args())
    main(args)

