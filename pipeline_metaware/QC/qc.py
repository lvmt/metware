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



import os
import sys
import textwrap
from collections import defaultdict

try:
    from Utils import utils
except Exception as e:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
    from Utils import utils



class Quality_Control:

    def __init__(self, args):
        self.projdir = args['projdir']
        self.sample_file = args['sample_file']
        self.analysis_list = args['analysis_list']
        self.fq_info = utils.get_fq_info(self.sample_file)


    def fastp(self, sampleID, fq1, fq2):
        # 适用于无宿主污染的环境样本
        cmd = textwrap.dedent(f'''
        /share/software/apps/anaconda2/bin/fastp \\
            -i {fq1} \\
            -o {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            -I {fq2} \\
            -O {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            -j {sampleID}.fastp.json &&\\
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc.sh'
        utils.write_cmd(cmd, shellname)


    def trim_and_bowtie2(self, sampleID, fq1, fq2):
        # 质控加上去除宿主
        cmd = textwrap.dedent(f'''
        # # trim质控
        # java -jar /lustrefs/share3/Bioinfo/lvmengting/pipeline/metagenomics/Trimmomatic-0.36/trimmomatic-0.36.jar \\
        #     PE  \\
        #     -threads 4 \\
        #     {fq1} \\
        #     {fq2} \\
        #     {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.1.fq.gz \\
        #     {self.projdir}/1.Clean/{sampleID}/{sampleID}.se.1.fq.gz \\
        #     {self.projdir}/1.Clean/{sampleID}/{sampleID}.pe.2.fq.gz \\
        #     {self.projdir}/1.Clean/{sampleID}/{sampleID}.se.2.fq.gz \\
        #     ILLUMINACLIP:/lustrefs/share3/Bioinfo/lvmengting/pipeline/metagenomics/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:40:15 \\
        #     LEADING:2 \\
        #     TRAILING:2 SLIDINGWINDOW:4:2 \\
        #     MINLEN:25 &&
        
        # 去除宿主
        /share/software/apps/bowtie2/2.3.4/bowtie2 \\
            --sensitive -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 \\
            -I 200 -X 400 \\
            -p 6 \\
            -x ~/pipeline/metagenomics/database/fasta/human/human \\
            -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            --un-conc-gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.gz
        
        # 结果文件软连
        # ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.1.gz {self.projdir}/1.Clean/{sampleID}.final.1.gz
        # ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.2.gz {self.projdir}/1.Clean/{sampleID}.final.2.gz
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc1.sh'
        utils.write_cmd(cmd, shellname)


    def fastqc(self, sampleID):
        # 数据质量统计分析
        cmd = textwrap.dedent(f'''
        fastqc {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz -t 4
        fastqc {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz -t 4
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/fastqc.sh'
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
        # utils.mkdirs(f'{self.projdir}/1.Clean/multiqc')

        for sampleID in self.fq_info:
            fq1, fq2 = self.fq_info[sampleID]
            sample_dir = f'{self.projdir}/1.Clean/{sampleID}'
            utils.mkdirs(sample_dir)
            
            # 质控程序选择
            if '1.1' in analysis_list:
                self.fastp(sampleID, fq1, fq2)
            elif '1.2' in analysis_list:
                self.trim_and_bowtie2(sampleID, fq1, fq2)

            # self.fastqc(sampleID)  # 生成质控信息
            self.stat_info(sampleID)
            # self.multiqc() 



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str)
    # 质控分析参数
    parser.add_argument('--host', help='宿主来源')
    
    args = vars(parser.parse_args())
    Quality_Control(args).start()

   
