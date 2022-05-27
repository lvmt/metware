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
from GenePrediction.gene_prediction import Prediction
from TaxAnnotation.taxonomy_annotation import Taxonomy



from collections import defaultdict


class MainPipeline:
    
    def __init__(self, args):
        self.args = args
        # self.projdir = args['projdir']
        # self.sample_file = args['sample_file']
        # self.analysis_list = args['analysis_list']


    def start(self):
        Quality_Control(self.args).start()
        Meta_Assembly(self.args).start()
        Prediction(self.args).start()
        Taxonomy(self.args).start()





if __name__ == '__main__':
    analysis_info = {
        '1.1': '>>> QC: 无宿主',
        '1.2': '>>> QC: 有宿主',
        '2': '>>> Assembly: 组装',
        '3': '>>> GenePredict 基因预测',
        '4': '>>> TaxAnnotation 物种预测',
        '5': '>>> FunAnnotation 功能注释'
    }
    import argparse
    parser = argparse.ArgumentParser(description='metagenomics pipeline')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str)
    # 质控分析参数
    parser.add_argument('--host', help='宿主来源')
    # diamond 参数
    parser.add_argument('--threshold', help='evalue阈值', type=float, default=0.00001)
    args = vars(parser.parse_args())
    print(args)
    
    print('\033[1;32m分析步骤\033[0m')
    for step in args['analysis_list'].split(';'):
        print(analysis_info.get(step))

    MainPipeline(args).start()
