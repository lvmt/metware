#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-05 09:29:04
'''


from ctypes import util
from os import sep
import textwrap
from Utils import utils
import re 
import pandas as pd  


class Prediction:

    def __init__(self, args, fq_info):
        self.args = args 
        self.fq_info = fq_info
        self.projdir = self.args['projdir']


    def genemark(self, sampleID):
        # 基于单个样本进行基因预测 
        cmd = textwrap.dedent(f'''
        # ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/gmhmmp \\
        #     -a -d -f G \\
        #     -m ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod \\
        #     -A {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.fa \\
        #     -D {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.fa \\
        #     -o {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.contigs.fa.gff \\
        #     {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa

        ## 对预测结果名称进行修改
        sed '/>/s/gene/{sampleID}_gene/' {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.fa \\
            | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.fa

        ## 100nt过滤
        ~/.conda/envs/python2_lmt/bin/seqkit seq -m 100 \\
            {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.fa \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.100nt.fa
        
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/{sampleID}/prediction.sh'
        utils.write_cmd(cmd, shellname)


    def prodigal(self, sampleID):
        # 基于单个样本的基因预测 
        pass 


    def gene_cluster(self):
        # 对预测结果进行聚类
        nucleotide_list = []
        protein_list = []
        for sampleID in self.fq_info:
            nucleotide_tmp  = f'{self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.fa'
            protein_tmp = f'{self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.fa'
            nucleotide_list.append(nucleotide_tmp)
            protein_list.append(protein_tmp)

        nucleotide = ' '.join(nucleotide_list)
        protein = ' '.join(protein_list)
        cmd = textwrap.dedent(f'''
        cat {nucleotide} > {self.projdir}/3.GenePrediction/Cluster/total.nucleotide.fa

        cat {protein} > {self.projdir}/3.GenePrediction/Cluster/total.protein.fa

        ~/.conda/envs/python2_lmt/bin/cd-hit \\
            -c 0.95 \\
            -aS 0.9 \\
            -G 0 \\
            -g 1 \\
            -d 0 \\
            -M 6000 \\
            -i {self.projdir}/3.GenePrediction/Cluster/total.nucleotide.fa \\
            -o {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa 

        # 翻译成蛋白序列
        seqkit translate \\
            {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            > {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.fa 

        bowtie2-build -f {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/Cluster/gene_cluster.1.sh'
        utils.write_cmd(cmd, shellname)
        

    def stat_abundance(self, sampleID):
        ref = f'{self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide'
        cmd = textwrap.dedent(f'''
        # bowtie2 比对
        bowtie2 -p 4 \\
            -x {ref} \\
            -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            --end-to-end --sensitive -I 200 -X 400 \\
            -S {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.sam  

        # gene支持数
        grep -v '^@'  {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.sam  \\
            | cut -f3 | awk '{{if($1 != "*")print $0}}'  \\
            | sort | uniq -c |  awk 'BEGIN{{FS=" "; OFS=","}}{{print $2, $1}}' \\
            | awk 'BEGIN{{FS=",";OFS=","}}{{if($2>2)print $1, $2;else print $1,"0"}}' \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.count.csv 

        # 统计基因长度信息
        bioawk -c fastx  '{{print $name, length($seq)}}' \\
            {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            | tr '\t' ',' > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.length.csv 

        python3 get_abundance_table.py \\
            --gene_count {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.count.csv  \\
            --gene_length {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.length.csv  \\
            --abundance_table {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.table
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.stat_abundance.sh'
        utils.write_cmd(cmd, shellname)


    def merge_abundance_table(self):
        # 合并全部样本的丰度表格
        table_list = []
        for sample in self.fq_info:
            table = f'{self.projdir}/3.GenePrediction/{sample}/{sample}.abundance.table'
            table_list.append(table)
        table_files = " ".join(table_list)
        cmd = textwrap.dedent(f'''
        python3 merge_abundance_table.py \\
            --abundance_tables {table_files} \\
            --merge_abundance  {self.projdir}/3.GenePrediction/Cluster/total.abundance.tables
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/Cluster/merge_all_abundance_tabls.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        utils.mkdirs(f'{self.projdir}/3.GenePrediction/Cluster')
        # 基因预测
        for sampleID in self.fq_info:
            utils.mkdirs(f'{self.projdir}/3.GenePrediction/{sampleID}')
            self.genemark(sampleID)

        # 聚类
        self.gene_cluster()
        
        # 丰度计算
        for sampleID in self.fq_info:
            self.stat_abundance(sampleID)
            
        self.merge_abundance_table()

        
    