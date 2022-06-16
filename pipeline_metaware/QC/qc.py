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


import yaml
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
        self.args = args 
        self.projdir = self.args['projdir']
        self.sample_file = self.args['sample_file']
        self.analysis_list = self.args['analysis_list']
        self.fq_info = utils.get_fq_info(self.sample_file)
        

    def link_rawdata(self):
        # 先软连原始数据
        for sample in self.fq_info:
            utils.mkdirs(f'{self.projdir}/0.RawData/{sample}')
            fq1, fq2 = self.fq_info[sample]
            os.system(f'ln -sf {fq1} {self.projdir}/0.RawData/{sample}/{sample}.fq.1.gz')
            os.system(f'ln -sf {fq2} {self.projdir}/0.RawData/{sample}/{sample}.fq.2.gz')


    def fastp(self, sampleID, fq1, fq2):
        # 适用于无宿主污染的环境样本
        cmd = textwrap.dedent(f'''
        /share/software/apps/anaconda2/bin/fastp \\
            -i {fq1} \\
            -o {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.1.gz \\
            -I {fq2} \\
            -O {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.2.gz \\
            -j {sampleID}.fastp.json &&\\
        
        # 结果文件软连
        ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.1.gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz
        ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.2.gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc.sh'
        utils.write_cmd(cmd, shellname)


    def qc_and_bowtie2(self, sampleID, fq1, fq2):
        # 质控加上去除宿主
        host = self.args['HOST'][self.args['host']]
        cmd = textwrap.dedent(f'''
        # 去除宿主
        /share/software/apps/bowtie2/2.3.4/bowtie2 \\
            --sensitive -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 \\
            -I 200 -X 400 \\
            -p 6 \\
            -x {host} \\
            -1 {fq1} \\
            -2 {fq2} \\
            --un-conc-gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.gz
                
        # fastqc质控
        /share/software/apps/anaconda2/bin/fastp \\
            -i {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.1.gz \\
            -o {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.1.gz \\
            -I {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.2.gz \\
            -O {self.projdir}/1.Clean/{sampleID}/{sampleID}.clean.fq.2.gz \\
            -j {sampleID}.fastp.json &&\\

        # 结果文件软连
        ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.1.gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz
        ln -sf {self.projdir}/1.Clean/{sampleID}/{sampleID}.remove_host.fq.2.gz {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/qc.sh'
        utils.write_cmd(cmd, shellname)


    def md5(self, sampleID):
        # 计算质控后clean fq的md5值
        cmd = textwrap.dedent(f'''
        md5sum {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz > md5.txt
        md5sum {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz >> md5.txt
        ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/md5.sh'
        utils.write_cmd(cmd, shellname)


    def stat_info(self, sampleID, fq1, fq2):
        '''
        # 条数/碱基数/Q20碱基数/Q30碱基数/GC碱基数
        '''
        cmd = textwrap.dedent(f'''
        # # 统计rawdata 和 clean data, 进行比对
        # # rawdata
        ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/fastq_stat \\
            {fq1} \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.raw.1.info.xls 

        # ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/fastq_stat \\
            {fq2} \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.raw.2.info.xls 
        
        # ## clean
        ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/fastq_stat \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.info.xls

        # ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/fastq_stat \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.info.xls

        ## 绘图
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/plot.py \\
            {self.projdir}/1.Clean/{sampleID}/{sampleID}.fastp.json
        
        # ## 质控信息统计
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/stat.py \\
            --raw1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.raw.1.info.xls  \\
            --raw2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.raw.2.info.xls  \\
            --clean1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.info.xls\\
            --clean2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.info.xls\\
            --stat_result {self.projdir}/1.Clean/{sampleID}/{sampleID}.stat.txt 

         ''')
        shellname = f'{self.projdir}/1.Clean/{sampleID}/stat.sh'
        utils.write_cmd(cmd, shellname)


    def merge_stat(self):
        stat_files = []
        for sample in self.fq_info:
            stat = f'{self.projdir}/1.Clean/{sample}/{sample}.stat.txt'
            stat_files.append(stat)
        stat_files = ' '.join(stat_files)
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/QC/bin/combine_stat.py \\
            --stat_files {stat_files}\\
            --stat_merge {self.projdir}/1.Clean/all.samples.stat.txt
        ''')
        shellname = f'{self.projdir}/1.Clean/merge_stat.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        self.link_rawdata()

        for sampleID in self.fq_info:
            fq1, fq2 = self.fq_info[sampleID]
            sample_dir = f'{self.projdir}/1.Clean/{sampleID}'
            utils.mkdirs(sample_dir)
            
            # 质控程序选择
            if '1.1' in self.analysis_list.split(','):
                self.fastp(sampleID, fq1, fq2)
            elif '1.2' in self.analysis_list.split(','):
                self.qc_and_bowtie2(sampleID, fq1, fq2)

            self.stat_info(sampleID, fq1, fq2)
            self.md5(sampleID)

        self.merge_stat()
        



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str)
    parser.add_argument('--host', help='宿主来源')
    parser.add_argument('--config', help='一些配置参数, 命令行参数优先级高于配置文件')
    
    args = vars(parser.parse_args())
    config_info = yaml.safe_load(open(args['config'])) 

    for key, value in args.items():
        if value:
            config_info.update({key: value})

    print(config_info)

    Quality_Control(config_info).start()
