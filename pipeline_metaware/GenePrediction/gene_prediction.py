#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-17 11:02:26
'''


import re 
import os
import sys
import textwrap
import pandas as pd  

try:
    from Utils import utils
except Exception as e:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
    from Utils import utils



class Prediction:

    def __init__(self, args):
        self.args = args 
        self.sample_file = self.args['sample_file']
        self.projdir = self.args['projdir']
        self.analysis_list = args['analysis_list']
        self.fq_info = utils.get_fq_info(self.sample_file)


    def prediction_unmap_reads(self):
        # 混合组装contig进行基因预测
        utils.mkdirs(f'{self.projdir}/3.GenePrediction/unmap_predict')
        cmd = textwrap.dedent(f'''
        ## 基因预测
        # ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/gmhmmp \\
        #     -a -d -f G \\
        #     -m ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod \\
        #     -A {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.fa \\
        #     -D {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.fa \\
        #     -o {self.projdir}/3.GenePrediction/unmap_predict/unmap.contigs.fa.gff \\
        #     {self.projdir}/2.Assembly/unmap_assembly/unmap_reads.contigs.fa &&\\

        # ## 对预测结果名称进行修改及100nt过滤: 核酸序列
        # sed '/>/s/gene/unmap_gene/' {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.fa \\
        #     | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \\
        #     > {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.rename.fa  &&\\

        # ~/.conda/envs/python2_lmt/bin/seqkit seq -m 100 \\
        #     {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.rename.fa \\
        #     > {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.rename.100nt.fa &&\\

        # ln -sf {self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.rename.100nt.fa \\
        #     {self.projdir}/3.GenePrediction/unmap_predict/unmap.cds.fa &&\\

        ## 获取长度大于100nt的基因名称
         ~/.conda/envs/python2_lmt/bin/seqkit seq -n \\
            {self.projdir}/3.GenePrediction/unmap_predict/unmap.cds.fa  \\
            > {self.projdir}/3.GenePrediction/unmap_predict/unmap.gene_name.100nt

        ## 对预测结果名称进行修改及100nt过滤: 蛋白序列
        sed '/>/s/gene/unmap_gene/' {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.fa \\
            | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \\
            > {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.rename.fa  &&\\

        ~/.conda/envs/python2_lmt/bin/seqkit  \\
            grep \\
            {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.rename.fa \\
            -f {self.projdir}/3.GenePrediction/unmap_predict/unmap.gene_name.100nt \\
            > {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.rename.100nt.fa &&\\

        ln -sf {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.rename.100nt.fa \\
            {self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.final.fa &&\\

        ## 对预测结果进行统计分析
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/stat_cds_info.py \\
            --fa {self.projdir}/3.GenePrediction/unmap_predict/unmap.cds.fa \\
            --result  {self.projdir}/3.GenePrediction/unmap_predict/unmap.cds.stat_fa  &&\\

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/summary_cds.py \\
            --infile {self.projdir}/3.GenePrediction/unmap_predict/unmap.cds.stat_fa  \\
            --result_suffix {self.projdir}/3.GenePrediction/unmap_predict/unmap
        
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/unmap_predict/prediction.sh'
        utils.write_cmd(cmd, shellname)


    def prediction_sample(self, sampleID):
        utils.mkdirs(f'{self.projdir}/3.GenePrediction/{sampleID}')
        # 基于单个样本进行基因预测 
        cmd = textwrap.dedent(f'''
        ## 基因预测
        # ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/gmhmmp \\
        #     -a -d -f G \\
        #     -m ~/pipeline/metagenomics/software/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod \\
        #     -A {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.fa \\
        #     -D {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.fa \\
        #     -o {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.contigs.fa.gff \\
        #     {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa  &&\\

        ## 对预测结果名称进行修改及100nt过滤: 核酸序列
        # sed '/>/s/gene/{sampleID}_gene/' {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.fa \\
        #     | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \\
        #     > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.fa  &&\\

        # ~/.conda/envs/python2_lmt/bin/seqkit seq -m 100 \\
        #     {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.fa \\
        #     > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.100nt.fa &&\\

        # ln -sf {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.100nt.fa \\
        #     {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.cds.fa &&\\

        ## 获取长度大于100nt的基因名称
        ~/.conda/envs/python2_lmt/bin/seqkit seq -n \\
            {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.cds.fa \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene_name.100nt &&\\
 
        ## 对预测结果名称进行修改及100nt过滤: 蛋白序列
        sed '/>/s/gene/{sampleID}_gene/' {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.fa \\
            | sed '/>/s/GeneMark.hmm//'| awk  -F '|' '{{print $1}}' \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.rename.fa  &&\\

        ~/.conda/envs/python2_lmt/bin/seqkit  \\
            grep \\
            {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.rename.fa \\
            -f {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene_name.100nt \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.rename.100nt.fa &&\\
        
        ln -sf {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.rename.100nt.fa \\
            {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.final.fa &&\\
        
        # 对预测结果进行统计分析
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/stat_cds_info.py \\
            --fa {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.cds.fa \\
            --result  {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.cds.stat_fa  &&\\

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/summary_cds.py \\
            --infile {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.cds.stat_fa \\
            --result_suffix {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}

        ''')
        shellname = f'{self.projdir}/3.GenePrediction/{sampleID}/prediction.sh'
        utils.write_cmd(cmd, shellname)


    def gene_cluster(self):
        # 对预测结果进行聚类
        utils.mkdirs(f'{self.projdir}/3.GenePrediction/Cluster')
        nucleotide_list = []
        protein_list = []
        for sampleID in self.fq_info:
            nucleotide_tmp  = f'{self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.nucleotide.rename.100nt.fa'
            protein_tmp = f'{self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.protein.rename.100nt.fa'
            nucleotide_list.append(nucleotide_tmp)
            protein_list.append(protein_tmp)

        # 加入unmap的预测结果
        nucleotide_list.append(f'{self.projdir}/3.GenePrediction/unmap_predict/unmap.nucleotide.rename.100nt.fa')
        protein_list.append(f'{self.projdir}/3.GenePrediction/unmap_predict/unmap.protein.rename.100nt.fa')
        nucleotide_files = ' '.join(nucleotide_list)
        protein_files = ' '.join(protein_list) 
        protein_tools = '/lustrefs/share3/Bioinfo/lvmengting/.conda/envs/python2_lmt/bin/cd-hit'
        nucleotide_tools = '/lustrefs/share3/Bioinfo/lvmengting/.conda/envs/python2_lmt/bin/cd-hit-est'
        files = {'protein': protein_files, 'nucleotide': nucleotide_files}
        tools = {'protein': protein_tools, 'nucleotide': nucleotide_tools}

        for _type in ['protein', 'nucleotide']:
            infiles = files[_type]
            tool = tools[_type]

            cmd = textwrap.dedent(f'''
            # cat {infiles} > {self.projdir}/3.GenePrediction/Cluster/total.{_type}.fa

            {tool} \\
                -c 0.95 \\
                -aS 0.9 \\
                -G 0 \\
                -g 1 \\
                -d 0 \\
                -M 100000 \\
                -i {self.projdir}/3.GenePrediction/Cluster/total.{_type}.fa \\
                -o {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.{_type}.fa &&\\

            # # 翻译成蛋白序列
            # seqkit translate \\
            #     {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            #     > {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.fa  &&\\

            # bowtie2-build -f {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            #     {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide
            ''')
            shellname = f'{self.projdir}/3.GenePrediction/Cluster/{_type}_cluster.sh'
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
            -S {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.sam  &&\\

        # gene支持数
        grep -v '^@'  {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.sam  \\
            | cut -f3 | awk '{{if($1 != "*")print $0}}'  \\
            | sort | uniq -c |  awk 'BEGIN{{FS=" "; OFS=","}}{{print $2, $1}}' \\
            | awk 'BEGIN{{FS=",";OFS=","}}{{if($2>2)print $1, $2;else print $1,"0"}}' \\
            > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.count.csv &&\\

        # 统计基因长度信息
        bioawk -c fastx  '{{print $name, length($seq)}}' \\
            {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            | tr '\t' ',' > {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.length.csv  &&\\

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/get_abundance_table.py \\
            --gene_count {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.count.csv  \\
            --gene_length {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.gene.length.csv  \\
            --abundance_table {self.projdir}/3.GenePrediction/{sampleID}/{sampleID}.abundance.table  &&\\
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
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/GenePrediction/bin/merge_abundance_table.py \\
            --abundance_tables {table_files} \\
            --merge_abundance  {self.projdir}/3.GenePrediction/Cluster/total.abundance.tables
        ''')
        shellname = f'{self.projdir}/3.GenePrediction/Cluster/merge_all_abundance_tabls.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        analysis_list = self.analysis_list.split(';')
        if '3' in analysis_list:
            # 单样本基因预测
            for sampleID in self.fq_info:                
                self.prediction_sample(sampleID)

            # unmap read 基因预测
            self.prediction_unmap_reads()

            # 聚类
            self.gene_cluster()
            
            # 丰度计算
            for sampleID in self.fq_info:
                self.stat_abundance(sampleID)
                
            self.merge_abundance_table()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str, default='3')
    
    args = vars(parser.parse_args())
    print(args)
    Prediction(args).start()

        