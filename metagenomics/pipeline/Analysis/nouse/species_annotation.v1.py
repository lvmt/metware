#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-09 15:58:23
'''

# 物种注释：基于diamond
'''
1 diamond比对
2 从输出结果中提取accession.version号
3 从../prot.accession2taxid.gz提取taxid
4 基于taxid,注释物种信息 
5 给通过accession.version映射基因与物种关系
'''


import textwrap
from Utils import utils


class LineageAnnotation:

    def __init__(self, args):
        self.args = args
        self.projdir = self.args['projdir']
        self.threshold = self.args['threshold']  # 默认值0.001


    def diamond(self):
        representive_fa = f'{self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa'
        blast_result = f'{self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.daa'
        cmd = textwrap.dedent(f'''
        # 物种注释
        ~/pipeline/metagenomics/software/diamond \\
            blastx \\
            -d bacteria.nr.dmnd \\
            -q {representive_fa } \\
            -o {blast_result} \\
            --outfmt 100 \\
            --evalue {self.threshold} &&\\
              
        # 获取accession.version
        cut -f2 {blast_result } | sort | uniq > {blast_result }.acc

        # 获取acc ~ taxid 对应关系
        zcat ~/pipeline/metagenomics/database/prot.accession2taxid.gz \\
            | csvtk -t grep -f 2 -P {blast_result }.acc \\
            | cut -f2,3 \\
            > {blast_result}.acc.taxid

        # LCA 分析, 每条序列返回唯一比对结果
        python3 lca.py \\
            --alignment {blast_result} \\
            --relation {blast_result}.acc.taxid \\
            --result {blast_result}.only_taxid

        # 根据TaxID获取完整谱系(lineage)
        sed '/accession/d' {blast_result}.only_taxid \\
            |cut -f2  \\
            | ~/.conda/envs/python3_lmt/bin/taxonkit lineage \\
            | taxonkit reformat -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F -P \\
            | csvtk cut -t -f -2 \\
            | csvtk add-header -t -n taxid,kindom,phylum,class,order,family,genus,species \\
            > {blast_result }.taxonomy

        # 获取taxonomy table, 将物种信息映射到序列上
        python3 get_taxonomy_table.py  \\
            --gene_taxid {blast_result}.only_taxid \\
            --taxid_taxonomy {blast_result }.taxonomy \\
            --taxonomy_anno_table {blast_result }.taxonomy_table 

        # 合并丰度和物种注释表格
        python3 combine_abundance_taxonomy.py \\
            --abundance \\
            --taxonomy \\
            --combine_table \\
        ''')
        shellname = f'{self.projdir}/4.GeneAnnotation/anno/diamond.sh'
        utils.write_cmd(cmd, shellname)


    def combine_gene_abundance(self):
        # 整合定量表和物种注释表
        gene_table = f'{self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.m8.taxonomy'
        abundance_table = f'{self.projdir}/3.GenePrediction/Cluster/total.abundance.tables'




    def start(self):
        self.diamond()
        
