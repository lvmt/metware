#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-17 08:51:13
'''

# 设置任务怎么跑
import os 
import subprocess


class Run:
    
    def __init__(self, args):
        self.args = args 
        self.projdir = self.args['projdir']
        self.analysis_list = self.args['analysis_list']
        self.sample_file = self.args['sample_file']
        self.samples = self.get_samples()
        

    
    def get_samples(self):
        samples = []
        with open(self.sample_file, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                samples.append(linelist[0])
        return samples
        

    def get_depend_ids(self, _ids):
        return ':'.join(_ids) if isinstance(_ids, list) else _ids


    def qc(self):
        print('>>> 任务投递 QC ')
        qc_jobids = []
        for sample in self.samples:
            qc = f'{self.projdir}/1.Clean/{sample}/qc.sh'
            md5 = f'{self.projdir}/1.Clean/{sample}/md5.sh'
            stat = f'{self.projdir}/1.Clean/{sample}/stat.sh'
            
            qc_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {qc} --mem 4gb')
            os.system(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {md5} --mem 2gb -W {qc_jobid}')
            stat_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {stat} --mem 4gb -W {qc_jobid}')
            qc_jobids.append(stat_jobid)

        ## 合并统计结果
        qc_stat_sh = f'{self.projdir}/1.Clean/merge_stat.sh'
        qc_jobids = ':'.join(qc_jobids)
        _id = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {qc_stat_sh} --mem 4gb -W {qc_jobids}')
        return _id 


    def assembly(self, qc_jobids):
        print('>>> 任务投递：Assembly')
        depend_ids = self.get_depend_ids(qc_jobids)
        sample_assembly_jobids = []

        for sample in self.samples:
            sample_assembly_sh = f'{self.projdir}/2.Assembly/{sample}/assembly.sh'
            if qc_jobids:
                jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {sample_assembly_sh} --mem 40gb --thread 10 -W {depend_ids}')
            else:
                jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {sample_assembly_sh} --mem 40gb --thread 10')
            sample_assembly_jobids.append(jobid)
        
        ## 混合组装, 依赖单样本组装
        unmap_assembly_sh = f'{self.projdir}/2.Assembly/unmap_assembly/unmap_assembly.sh'
        unmap_depend = self.get_depend_ids(sample_assembly_jobids)
        unmap_assembly_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {unmap_assembly_sh} --mem 80gb --thread 20 -W {unmap_depend}')

        # 评估结果整理 
        merge_stat_sh = f'{self.projdir}/2.Assembly/merge_stat.sh'
        _id = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {merge_stat_sh} --mem 1gb -W {unmap_assembly_jobid}')
        return _id 

    
    def predict(self, assembly_jobids):
        pass 


    def taxon(self, predict_jodids):
        pass 

    def func(self, taxon_jobids):
        pass 


    def start(self):
        qc_jobids = ''
        assembly_jobids = ''
        predict_jobids = ''
        taxo_jobids = ''
        func_jboids = ''

        analysis_list = [ana.split('.')[0] for ana in self.analysis_list.split(',')]
        if '1' in analysis_list:
            qc_jobids = self.qc()

        if '2' in analysis_list:
            assembly_jobids = self.assembly(qc_jobids)
