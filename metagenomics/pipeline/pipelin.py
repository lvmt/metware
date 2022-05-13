#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-26 14:26:02
@Last Modified by:   lvmengting
@Last Modified time: 2022-04-26 14:26:02
'''


from Utils import utils
from QC.qc import Quality_Control
from Assembly.assembly import Meta_Assembly
from Analysis.gene_prediction import Prediction
from Analysis.species_annotation import LineageAnnotation


from collections import defaultdict


class MainPipeline:
    
    def __init__(self, args):
        self.args = args
        self.projdir = args['projdir']
        self.sample_file = args['sample_file']
        self.analysis_list = args['analysis_list']


    def get_fq_info(self):
        '''
        sample  type fq1    fq2
        '''
        fq_info_dict = defaultdict(list)
        with open(self.sample_file, 'r') as fr:
            for line in fr:
                if not line.strip(''):
                    continue
                linelist = line.strip('\n').split('\t')
                if linelist[0].lower() in ('sample'):
                    continue
                fq_info_dict[linelist[0]] = [linelist[2], linelist[3]]
        
        return fq_info_dict


    def start(self):
        fq_info = self.get_fq_info()

        #Quality_Control(self.args, fq_info).start()
        # Meta_Assembly(args, fq_info).start()
        # Prediction(args, fq_info).start()
        LineageAnnotation(args).start()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str)

    # diamond 参数
    parser.add_argument('--threshold', help='evalue阈值', type=float, default=0.001)

    args = vars(parser.parse_args())

    MainPipeline(args).start()

    print(args)






