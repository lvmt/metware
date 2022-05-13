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
from random import sample
from Utils import utils
import textwrap



class Meta_Assembly:
    
    def __init__(self, args, fq_info):
        self.args = args
        self.projdir = args['projdir']
        self.analysis_list = args['analysis_list']
        self.fq_info = fq_info


    def megahit(self, sampleID):
        # 占用内存小速度快,组装结果略微粗糙
        cmd = textwrap.dedent(f'''
        # ~/.conda/envs/python2_lmt/bin/megahit \\
        # -t 4 \\
        # -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
        # -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
        # --out-dir {self.projdir}/2.Assembly/{sampleID}/megahit \\
        # --out-prefix {sampleID} \\
        # --k-min 35 \\
        # --k-max 95 \\
        # --k-step 20 \\
        # --min-contig-len 500 \\
        # -f 

        ln -sf {self.projdir}/2.Assembly/{sampleID}/megahit/{sampleID}.contigs.fa \\
            {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa

        # 获取未被利用的reads
        # /share/software/apps/anaconda2/bin/bowtie2建立所有存在问题
        /share/software/apps/bowtie2/2.3.4/bowtie2-build \\
            -f {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa \\
            {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs \\
            --threads 2 --seed 1
        
        bowtie2 -I 200 -X 400 \\
            -x {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs \\
            -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
            -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
            --un-conc-gz {self.projdir}/2.Assembly/{sampleID}/{sampleID}.unmaped.fq.gz
        ''')

        shellname = f'{self.projdir}/2.Assembly/{sampleID}/assembly.sh'
        utils.write_cmd(cmd, shellname)


    def metaspades(self, sampleID):
        # 精细拼接
        cmd =  textwrap.dedent(f'''
        ~/pipeline/metagenomics/software/SPAdes-3.13.0-Linux/bin/metaspades.py \\
        -t 4 -m 20 \\
        -1 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.1.gz \\
        -2 {self.projdir}/1.Clean/{sampleID}/{sampleID}.fq.2.gz \\
        -o {self.projdir}/2.Assembly/{sampleID}
        ''')
        shellname = f'{self.projdir}/2.Assembly/{sampleID}/{sampleID}.assembly.sh'
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
        ~/.conda/envs/python2_lmt/bin/megahit \\
        -t 6 \\
        -1 {unmap1} \\
        -2 {unmap2} \\
        --out-dir {self.projdir}/2.Assembly/unmap_assembly/megahit \\
        --out-prefix unmap_reads \\
        --k-min 35 \\
        --k-max 95 \\
        --k-step 20 \\
        --min-contig-len 500 \\
        -f 

        ln -sf {self.projdir}/2.Assembly/unmap_assembly/megahit/unmap_reads.contigs.fa \\
            {self.projdir}/2.Assembly/unmap_assembly/unmap_reads.contigs.fa
        
        ''')
        shellname = f'{self.projdir}/2.Assembly/unmap_assembly/unmap_assembly.sh'
        utils.write_cmd(cmd, shellname)


    def assembly_assessment(self, sampleID):
        # 对组装的结果进行评估
        cmd = textwrap.dedent(f'''
        python3 ~/pipeline/metagenomics/software/quast-4.6.3/quast.py \\
            {self.projdir}/2.Assembly/{sampleID}/{sampleID}.contigs.fa \\
            -o {self.projdir}/2.Assembly/{sampleID}/quast \\
            -t 4 
        ''')
        shellname = f'{self.projdir}/2.Assembly/{sampleID}/assessment.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        analysis_list = self.analysis_list.split(';')

        for sampleID in self.fq_info:
            utils.mkdirs(f'{self.projdir}/2.Assembly/{sampleID}')
            if '2.1' in analysis_list:
                self.megahit(sampleID)
            elif '2.2' in analysis_list:
                self.metaspades(sampleID)

            self.assembly_assessment(sampleID)
            self.assembly_unmap_reads()


