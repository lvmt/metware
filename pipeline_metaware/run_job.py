#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-17 08:51:13
'''

# 设置任务怎么跑
import os 
import subprocess
import glob 


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


    def qc(self):  # sourcery skip: inline-immediately-returned-variable
        print('>>> 任务投递 QC ')
        qc_jobids = []
        for sample in self.samples:
            qc = f'{self.projdir}/1.Clean/{sample}/qc.sh'
            md5 = f'{self.projdir}/1.Clean/{sample}/md5.sh'
            stat = f'{self.projdir}/1.Clean/{sample}/stat.sh'
            
            qc_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {qc} --mem 4gb')
            subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {md5} --mem 2gb -W {qc_jobid}')
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
            jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {sample_assembly_sh} --mem 40gb --threads 10 -W {depend_ids}')
            sample_assembly_jobids.append(jobid)

        ## 混合组装, 依赖单样本组装
        unmap_assembly_sh = f'{self.projdir}/2.Assembly/unmap_assembly/unmap_assembly.sh'
        unmap_depend = self.get_depend_ids(sample_assembly_jobids)
        unmap_assembly_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {unmap_assembly_sh} --mem 80gb --threads 20 -W {unmap_depend}')

        # 评估结果整理 
        merge_stat_sh = f'{self.projdir}/2.Assembly/merge_stat.sh'
        return subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {merge_stat_sh} --mem 1gb -W {unmap_assembly_jobid}') 

    
    def predict(self, assembly_jobids):
        # sourcery skip: inline-immediately-returned-variable, merge-list-append, move-assign-in-block
        '''
        1. 先基因预测
        2. 聚类去冗余
        3. 丰度统计
        '''
        print('>>> 任务投递：Predict')
        depend_ids = self.get_depend_ids(assembly_jobids)
        ## 基因预测
        predict_jobids = []
        unmap_predict_sh = f'{self.projdir}/3.GenePrediction/unmap_predict/prediction.sh'
        jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {unmap_predict_sh} --mem 80gb --threads 20 -W {depend_ids}')
        predict_jobids.append(jobid)

        for sample in self.samples:
            predict_sh = f'{self.projdir}/3.GenePrediction/{sample}/prediction.sh'
            jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {predict_sh} --mem 40gb --threads 10 -W {depend_ids}')
            predict_jobids.append(jobid)

        ## 聚类 
        cluster_depend_ids = self.get_depend_ids(predict_jobids)
        cluster_sh = f'{self.projdir}/3.GenePrediction/Cluster/protein_cluster.sh'
        cluster_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {cluster_sh} --mem 100gb -q middle2 --threads 40 -W {cluster_depend_ids}')

        nucle_cluster_sh = f'{self.projdir}/3.GenePrediction/Cluster/nucleotide_cdhit.sh'
        nucle_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {nucle_cluster_sh} --mem 40gb --threads 10 -W {cluster_jobid}')

        ## 丰度统计
        abundance_jobids = []
        for sample in self.samples:
            stat_sh = f'{self.projdir}/3.GenePrediction/{sample}/stat_abundance.sh'
            jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {stat_sh} --mem 40gb  --threads 10 -W {nucle_jobid}')
            abundance_jobids.append(jobid)

        ## 丰度汇总统计
        merge_depend_ids = self.get_depend_ids(abundance_jobids)
        merge_stat_sh = f'{self.projdir}/3.GenePrediction/Cluster/merge_all_abundance_table.sh'
        merge_stat_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {merge_stat_sh} --mem 5gb -W {merge_depend_ids}')

        ## 绘图
        plot_sh = f'{self.projdir}/3.GenePrediction/Cluster/plot.sh'
        _id = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {plot_sh} --mem 5gb -W {merge_stat_jobid}')
        return _id 


    def taxon(self, predict_jodids):
        # sourcery skip: inline-immediately-returned-variable
        print('>>> 任务投递：Taxonomy')
        depend_ids = self.get_depend_ids(predict_jodids)
        # diamond比对
        diamond_sh = f'{self.projdir}/4.TaxAnnotation/diamond.sh'
        diamond_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {diamond_sh} --mem 100gb -q middle2 --threads 40 -W {depend_ids}')
        # megan
        megan_sh = f'{self.projdir}/4.TaxAnnotation/megan.sh'
        megan_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {megan_sh} --mem 40gb --threads 10 -W {diamond_jobid}')
        # merge taxonomy and abundance 
        merge_sh = f'{self.projdir}/4.TaxAnnotation/merge_abundance_taxonomy.sh'
        merge_jobid = subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {merge_sh} --mem 10gb -W {megan_jobid}')
        # 绘图
        plot_sh_list = glob.glob(f'{self.projdir}/4.TaxAnnotation/Plot/*sh')
        for plot_sh in plot_sh_list:
            subprocess.getoutput(f'/lustrefs/share3/Bioinfo/lvmengting/tools/qsub.py {plot_sh} --mem 10gb --threads 2 -W {merge_jobid}')
        return merge_jobid 


    def func(self, taxon_jobids):
        pass 


    def start(self):
        # sourcery skip: assign-if-exp, introduce-default-else, move-assign-in-block
        qc_jobids = 'none'
        assembly_jobids = 'none'
        predict_jobids = 'none'
        taxo_jobids = 'none'
        func_jboids = 'none'

        analysis_list = [ana.split('.')[0] for ana in self.analysis_list.split(',')]
        if '1' in analysis_list:
            qc_jobids = self.qc()

        if '2' in analysis_list:
            assembly_jobids = self.assembly(qc_jobids)

        if '3' in analysis_list:
            predict_jobids = self.predict(assembly_jobids)

        if '4' in analysis_list:
            taxon_jobids = self.taxon(predict_jobids)


