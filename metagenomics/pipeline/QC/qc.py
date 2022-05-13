#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-25 16:21:30
'''


'''
质控模块：
基于配置文件， 生成对应的样本分析目录及脚本
'''

from Utils import utils

import os
import textwrap
from collections import defaultdict



class Quality_Control:

    def __init__(self, args, fq_info_dict):
        self.projdir = args['projdir']
        self.analysis_list = args['analysis_list']
        self.fq_info = fq_info_dict
        

    def fastqc(self, sampleID):
        # 数据质量统计分析
        cmd = textwrap.dedent(f'''
        fastqc {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz -t 4
        fastqc {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz -t 4
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/fastqc.sh'
        utils.write_cmd(cmd, shellname)

    
    def fastp(self, sampleID, fq1, fq2):
        # 适用于无宿主污染的环境样本
        cmd = textwrap.dedent(f'''
        /share/software/apps/anaconda2/bin/fastp \\
            -i {fq1} \\
            -o {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            -I {fq2} \\
            -O {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            -z 4
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc.sh'
        utils.write_cmd(cmd, shellname)


    def kneaddata(self, sampleID, fq1, fq2):
        # 集质控与去宿主一体
        return textwrap.dedent(f'''
        ~/.conda/envs/knead_lmt/bin/kneaddata \\
            -i {fq1} -i {fq1} \\
            -o {self.projdir}/1.Clean/{sampleID} \\
            -v -t 8 --remove-intermediate-output \\
            --trimmomatic ~/.conda/envs/knead_lmt/bin/trimmomatic \\
            --trimmomatic-options "ILLUMINACLIP:/lustrefs/share3/Bioinfo/lvmengting/pipeline/metagenomics/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50" \\
            --reorder --bowtie2-options "--very-sensitive --dovetail" \\
            -db ~/pipeline/metagenomics/database/fasta/tair10
        ''')


    def trim_and_bowtie2(self, sampleID, fq1, fq2):
        # 质控加上去除宿主
        cmd = textwrap.dedent(f'''
        java -jar /lustrefs/share3/Bioinfo/lvmengting/pipeline/metagenomics/Trimmomatic-0.36/trimmomatic-0.36.jar \\
            PE  \\
            -threads 4 \\
            {fq1} \\
            {fq2} \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.1.fq.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.se.1.fq.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.2.fq.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.se.2.fq.gz \\
            ILLUMINACLIP:/lustrefs/share3/Bioinfo/lvmengting/pipeline/metagenomics/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:40:15 \\
            LEADING:2 \\
            TRAILING:2 SLIDINGWINDOW:4:2 \\
            MINLEN:25 &&

        bowtie2 \\
            --very-sensitive-local --no-unal -I 200 -X 400 -p 6 \\
            -x ~/pipeline/metagenomics/database/fasta/tair10 \\
            -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.1.fq.gz \\
            -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.2.fq.gz \\
            --un-conc-gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.gz 
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc.sh'
        utils.write_cmd(cmd, shellname)


    def multiqc(self):
        file_qc = f'{self.projdir}/1.Clean/multiqc/all_qc_samples.txt'
        with open(file_qc, 'w') as fw:
            for sampleID in self.fq_info:
                path = f'{self.projdir}/1.Clean/{sampleID}\n'
                fw.write(path)

        cmd = textwrap.dedent(f'''
        /lustrefs/share3/Bioinfo/lvmengting/.conda/envs/python3_lmt/bin/multiqc \\
            --file-list {file_qc}
        ''')
        shellname = f'{self.projdir}/1.Clean/multiqc/multiqc.sh'
        utils.write_cmd(cmd, shellname)


    def stat_info(self, sampleID):
        '''
        # 条数/碱基数/Q20碱基数/Q30碱基数/GC碱基数
        '''
        cmd = textwrap.dedent(f'''
        ~/gitlab/meta_genomics/metagenomics/metagenomics/pipeline/QC/fastq_stat \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.1.info.xls

        ~/gitlab/meta_genomics/metagenomics/metagenomics/pipeline/QC/fastq_stat \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.2.info.xls
         ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/stat.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        analysis_list = self.analysis_list.split(';')
        utils.mkdirs(f'{self.projdir}/1.Clean/multiqc')

        for sampleID in self.fq_info:
            fq1, fq2 = self.fq_info[sampleID]
            sample_dir = f'{self.projdir}/1.Clean/{sampleID}'
            utils.mkdirs(sample_dir)
            
            # 质控程序选择
            if '1.1' in analysis_list:
                self.fastp(sampleID, fq1, fq2)
            elif '1.2' in analysis_list:
                self.trim_and_bowtie2(sampleID, fq1, fq2)

            self.fastqc(sampleID)  # 生成质控信息
            self.stat_info(sampleID)
            self.multiqc() 

