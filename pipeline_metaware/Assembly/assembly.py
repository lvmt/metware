#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-26 14:16:55
@Last Modified by:   lvmengting
@Last Modified time: 2022-04-26 14:16:55
'''


'''
宏基因组组装步骤
'''

from ctypes import util
import sys 
import os 
import textwrap

from sklearn.cluster import MeanShift


try:
    from Utils import utils
except Exception as e:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
    from Utils import utils



class Meta_Assembly:
    
    def __init__(self, args):
        self.args = args
        self.projdir = args['projdir']
        self.analysis_list = args['analysis_list']
        self.sample_file = args['sample_file']
        self.fq_info = utils.get_fq_info(self.sample_file)


    def assemby_sample(self, sampleID):
        # 占用内存小速度快,组装结果略微粗糙
        utils.mkdirs(f'{self.projdir}/2.Assembly/{sampleID}')
        cmd = textwrap.dedent(f'''
        ## 单样本组装
        # ~/.conda/envs/python2_lmt/bin/megahit \\
        # -t 4 \\
        # -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz \\
        # -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz \\
        # --out-dir {self.projdir}/2.Assembly/{sampleID}/megahit \\
        # --out-prefix {sampleID} \\
        # --k-min 35 \\
        # --k-max 95 \\
        # --k-step 20 \\
        # --min-contig-len 500 \\
        # --keep-tmp-files \\
        # --continue  \\
        # -f &&\\

        # ln -sf {self.projdir}/2.Assembly/{sampleID}/megahit/{sampleID}.contigs.fa \\
        #     {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa &&\\

        ## 组装结果评估
        python3 ~/pipeline/metagenomics/software/quast-4.6.3/quast.py \\
            {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa \\
            -o {self.projdir}/2.Assembly/{sampleID}/quast \\
            --contig-thresholds 0,500,1000,5000,10000,25000,50000 \\
            --plots-format png

        ## 获取未被利用的reads, 后续参与混合组装
        # # /share/software/apps/anaconda2/bin/bowtie2建立所有存在问题
        # /share/software/apps/bowtie2/2.3.4/bowtie2-build \\
        #     -f {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa \\
        #     {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs \\
        #     --seed 1  &&\\
        
        /share/software/apps/bowtie2/2.3.4/bowtie2 -I 200 -X 400 \\
            -x {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs \\
            -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.1.gz \\
            -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.final.2.gz \\
            --un-conc-gz {self.projdir}/2.Assembly/{sampleID}/{sampleID}.unmaped.fq.gz
        ''')

        shellname = f'{self.projdir}/2.Assembly/{sampleID}/assembly.sh'
        utils.write_cmd(cmd, shellname)
        
    
    def assembly_unmap_reads(self):
        # 组装全部的未比对上的reads
        # 将结果输出单独的目录下unmap_assembly
        utils.mkdirs(f'{self.projdir}/2.Assembly/unmap_assembly')
        unmap1 = []
        unmap2 = []
        for sampleID in self.fq_info:
            tmp1 = f'{self.projdir}/2.Assembly/{sampleID}/{sampleID}.unmaped.fq.1.gz'
            tmp2 = f'{self.projdir}/2.Assembly/{sampleID}/{sampleID}.unmaped.fq.2.gz'
            unmap1.append(tmp1)
            unmap2.append(tmp2)

        unmap1 = ','.join(unmap1)
        unmap2 = ','.join(unmap2)
        
        cmd = textwrap.dedent(f'''
        # ## unmap reads混合组装
        # ~/.conda/envs/python2_lmt/bin/megahit \\
        # -t 6 \\
        # -1 {unmap1} \\
        # -2 {unmap2} \\
        # --out-dir {self.projdir}/2.Assembly/unmap_assembly/megahit \\
        # --out-prefix unmap_reads \\
        # --k-min 35 \\
        # --k-max 95 \\
        # --k-step 20 \\
        # --min-contig-len 500 \\
        # -f  &&\\

        # ln -sf {self.projdir}/2.Assembly/unmap_assembly/megahit/unmap_reads.contigs.fa \\
        #     {self.projdir}/2.Assembly/unmap_assembly/unmap_reads.contigs.fa
        
        ## 组装结果评估
        python3 ~/pipeline/metagenomics/software/quast-4.6.3/quast.py \\
            {self.projdir}/2.Assembly/unmap_assembly/unmap_reads.contigs.fa \\
            -o {self.projdir}/2.Assembly/unmap_assembly/quast \\
            --contig-thresholds 0,500,1000,5000,10000,25000,50000 \\
            --plots-format png
        ''')
        shellname = f'{self.projdir}/2.Assembly/unmap_assembly/unmap_assembly.sh'
        utils.write_cmd(cmd, shellname)


    def merge_stat(self):
        # 合并所有样本的统计结果
        stat_files = []
        for sample in self.fq_info:
            stat = f'{self.projdir}/2.Assembly/{sample}/quast/report.tsv'
            stat_files.append(stat)

        # 添加unmap reads的结果
        stat_files.append(f'{self.projdir}/2.Assembly/unmap_assembly/quast/report.tsv')
        stat_files = ' '.join(stat_files)

        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/Assembly/bin/merge_stat.py \\
            --stat_files {stat_files} \\
            --merge_result {self.projdir}/2.Assembly/all.samples.assembly_assessment.txt
        ''')
        shellname = f'{self.projdir}/2.Assembly/merge_stat.sh'
        utils.write_cmd(cmd, shellname)    


    def start(self):
        analysis_list = self.analysis_list.split(';')
        if '2' in analysis_list:
            for sampleID in self.fq_info:
                self.assemby_sample(sampleID)

            self.assembly_unmap_reads() # 混合组装
            self.merge_stat()   # 统计信息汇总




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str, default='2')
    
    args = vars(parser.parse_args())
    Meta_Assembly(args).start()
