#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-11 13:50:19
'''

'''
基于abundance.table和taxonomy中的gene字段,
获取每个样本每个物种的abundance内容
获取每个样本每个物种的基因数目
将上面的总表,拆成界门纲目科属种


物种:
    样本1:
        gene: [gene1, gene2]
        abundance: nun++
'''

import re
import pandas as pd  
import sys  
import pprint

from collections import defaultdict



def get_stat_info(abundance_table, taxonomy_table):
    df_abun = pd.read_csv(abundance_table, sep='\t')
    df_tax = pd.read_csv(taxonomy_table, sep='\t', header=None, names=['gene', 'taxonomy'])
    merge_df = pd.merge(df_abun, df_tax, on='gene', how='outer').fillna('k__other;p__other;c__other;o__other;f__other;g__other;s__other')
    merge_df.to_csv('gene_abundance_taxonomy.xls', index=None, sep='\t')
    columns = list(merge_df.columns)
    index_info = dict(zip(columns, range(len(columns))))
    del index_info['gene']
    del index_info['taxonomy']

    stat_info = defaultdict(dict)

    for row in merge_df.values:
        gene = row[0]
        tax = row[-1]
        if tax not in stat_info:
            stat_info[tax] = defaultdict(dict)
            for sample in index_info:
                stat_info[tax][sample]['gene'] = []
                stat_info[tax][sample]['abundance'] = 0
        else:
            for sample in index_info:
                num = float(row[index_info[sample]]) # 该样本的丰度数据 
                if num > 0:
                    stat_info[tax][sample]['gene'].append(gene)
                    stat_info[tax][sample]['abundance'] += num

    return stat_info, index_info


def get_single_tax_gene(samples, stat_info):
    ## 获取每个样本的物种丰度,基因数目,基因名称表格
    for sample in samples:
        with open(sample+'.taxonomy_abundance_geneNums.table', 'w') as fw:
            fw.write('taxonomy\tabundance\tgene_nums\tgene_names\n')
            for tax in stat_info:
                gene = stat_info[tax][sample]['gene']
                abundance = stat_info[tax][sample]['abundance']
                gene_nums = len(gene)
                gene_names = ';'.join(gene)
                fw.write(f'{tax}\t{abundance}\t{gene_nums}\t{gene_names}\n')


## 基因数目样本汇总表格, 丰度汇总表格
def get_combine_tax_gene_table(samples, stat_info):
    with open('taxonomy.geneNums.table', 'w') as fwg, open('taxonomy.abundance.table', 'w') as fwb:
        fwg.write('taxonomy\t{0}\n'.format('\t'.join(samples)))
        fwb.write('taxonomy\t{0}\n'.format('\t'.join(samples)))
        for tax in stat_info:
            tmp_g = [tax]
            tmp_b = [tax]
            for sample in samples:
                gene_nums = len(stat_info[tax][sample]['gene'])
                abundance = stat_info[tax][sample]['abundance']
                tmp_g.append(gene_nums)
                tmp_b.append(abundance)
            fwg.write('{0}\n'.format('\t'.join(map(str, tmp_g))))
            fwb.write('{0}\n'.format('\t'.join(map(str, tmp_b))))


def get_subclassify(samples, stat_info):
    # 拆分成界门纲目科属种
    '''
    k:
        tax:
            sample: nums / [gene]
    '''
    abundance_info = defaultdict(dict)  # 分组
    gene_info = defaultdict(dict)   # 分组
    total_abundance_info = defaultdict(dict)  # 物种总的丰度, 不分组
    total_gene_info = defaultdict(dict) # gene数目统计, 不分组 

    for _class in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
        abundance_info[_class] = defaultdict(dict)
        gene_info[_class] = defaultdict(dict)
        total_abundance_info[_class] = defaultdict(int)
        total_gene_info[_class] = defaultdict(list)

    for tax in stat_info:
        tax_list = tax.split(';')
        k = tax_list[0]
        p = ';'.join(tax_list[:2])
        c = ';'.join(tax_list[:3])
        o = ';'.join(tax_list[:4])
        f = ';'.join(tax_list[:5])
        g = ';'.join(tax_list[:6])
        s = ';'.join(tax_list[:7])
        relation = {'k': k, 'p': p, 'c': c, 'o': o, 'f': f, 'g': g, 's': s}

        # 获取分组结果
        for _class in relation:
            for sample in samples:
                abundance_num = stat_info[tax][sample]['abundance']
                gene_num = len(stat_info[tax][sample]['gene'])

                if relation[_class] not in abundance_info[_class]:
                    abundance_info[_class][relation[_class]] = defaultdict(int)
                    gene_info[_class][relation[_class]] = defaultdict(int)
                abundance_info[_class][relation[_class]][sample] += abundance_num
                gene_info[_class][relation[_class]][sample] += gene_num

        # 获取汇总结果
        for _class in relation:
            for sample in samples:
                total_abundance_info[_class][relation[_class]] += stat_info[tax][sample]['abundance']
                total_gene_info[_class][relation[_class]].extend(stat_info[tax][sample]['gene'])

    for _class in abundance_info:
        with open('taxonomy.abundance.table.' + _class, 'w') as fw:
            fw.write('class\t{0}\n'.format('\t'.join(samples)))
            for tax in abundance_info[_class]:
                fw.write(f'{tax}')
                for sample in samples:
                    nums = abundance_info[_class][tax][sample]
                    fw.write(f'\t{nums}')
                fw.write('\n')

    for _class in gene_info:
        with open('taxonomy.geneNums.table.' + _class, 'w') as fw:
            fw.write('class\t{0}\n'.format('\t'.join(samples)))
            for tax in gene_info[_class]:
                fw.write(f'{tax}')
                for sample in samples:
                    nums = gene_info[_class][tax][sample]
                    fw.write(f'\t{nums}')
                fw.write('\n')

    for _class in total_abundance_info:
        with open('total.taxonomy.abundance.table.' + _class, 'w') as fw:
            fw.write(f'class\tabundance\n')
            for tax in total_abundance_info[_class]:
                abundance = total_abundance_info[_class][tax]
                fw.write(f'{tax}\t{abundance}\n')


        with open('total.gene.geneNums.table.' + _class, 'w') as fw:
            fw.write('class\tnums\tgenes\n')
            for tax in total_gene_info[_class]:
                gene_len = len(set(total_gene_info[_class][tax]))
                genes = ';'.join(set(total_gene_info[_class][tax]))
                fw.write(f'{tax}\t{gene_len}\t{genes}\n')





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='表格统计')
    parser.add_argument('--abundance_table', help='输入文件,丰度表格')
    parser.add_argument('--taxonomy_table', help='物种注释表格')
    parser.add_argument('--out_dir', help='输出目录')

    args = vars(parser.parse_args())

    abun_table = args['abundance_table'] 
    tax_table = args['taxonomy_table'] 
    stat_info, index_info = get_stat_info(abun_table, tax_table)
    samples = index_info.keys()
    # pprint.pprint(stat_info)

    get_single_tax_gene(samples, stat_info)
    get_combine_tax_gene_table(samples, stat_info)
    get_subclassify(samples, stat_info)
