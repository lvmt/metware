#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-26 14:26:02
'''


import yaml
import pprint
from collections import defaultdict


from Utils import utils
from QC.qc import Quality_Control
from Assembly.assembly import Meta_Assembly
from GenePrediction.gene_prediction import Prediction
from TaxAnnotation.taxonomy_annotation import Taxonomy
from FunctionAnnotation.function_annotation import Function

from Result.release import Release
from run_job import Run


docs = """
\033[1;32m
                 _                                         _          
                | |                                       (_)         
  _ __ ___   ___| |_ __ _  __ _  ___ _ __   ___  _ __ ___  _  ___ ___ 
 | '_ ` _ \ / _ \ __/ _` |/ _` |/ _ \ '_ \ / _ \| '_ ` _ \| |/ __/ __|
 | | | | | |  __/ || (_| | (_| |  __/ | | | (_) | | | | | | | (__\__ \\
 |_| |_| |_|\___|\__\__,_|\__, |\___|_| |_|\___/|_| |_| |_|_|\___|___/
                           __/ |                                      
                          |___/                                       
                                \033[1;31m--pipeline for metagenomics analysis                                                                               
\033[0m
"""




class MainPipeline:
    
    def __init__(self, args):
        self.args = args
        self.config_info = yaml.safe_load(open(args['config'])) if args['config'] else {}
        self.update()
        # print(self.config_info)


    def update(self):
        # 将args的非空参数更新至config文件中
        for key, value in self.args.items():
            if value:
                self.config_info.update({key: value})
        return self


    def start(self):
        analysis_print = {
            '1.1': '>>> QC: 无宿主',
            '1.2': '>>> QC: 有宿主',
            '2': '>>> Assembly: 组装',
            '3': '>>> GenePredict 基因预测',
            '4': '>>> TaxAnnotation 物种预测',
            '5': '>>> FunAnnotation 功能注释'
        }

        # 定义运行逻辑
        analysis_info = {
            '1': Quality_Control,
            '2': Meta_Assembly,
            '3': Prediction,
            '4': Taxonomy,
            '5': Function
        }

        # 项目配置参数check
        pprint.pprint(self.config_info)

        # 只生成脚本，并不运行
        for ana in self.config_info['analysis_list'].split(','):
            print(f'\033[1;32m分析步骤 {analysis_print[ana]}\033[0m')
            ana = ana.split('.')[0]
            analysis_info[ana](self.config_info).start()

        # 运行脚本
        if self.config_info.get('run_job'):
            Run(self.config_info).start()
        
        
        




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=docs, formatter_class=argparse.RawTextHelpFormatter)
    ## 共同参数 
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str)
    parser.add_argument('--config', help='一些配置参数, 命令行的参数的优先级高于配置文件优先级')
    parser.add_argument('--run_job', action='store_true', default=False, help='是否投递任务,默认不投递')
    ## 质控分析参数
    parser.add_argument('--host', help='宿主来源, 用于去除对应的宿主')
    ##  物种注释参数 
    parser.add_argument('--nr_database', help='nr数据库路径')
    parser.add_argument('--nr_diamond_evalue', help='物种注释步骤,数据库比对阈值', default='0.00001')
    parser.add_argument('--mega_database', help='mega 数据库连接')
    parser.add_argument('--taxid_taxonomy_database', help='自建库, taxid对应的物种注释信息') 
    parser.add_argument('--metastat_compare', action='append', help='比较分组')
    parser.add_argument('--random_forest_compare', action='append', help='比较分组')
    ## 功能注释参数 
    parser.add_argument('--keeg_diamond_database', help='kegg蛋白序列数据库')
    parser.add_argument('--kegg_diamond_evalue', help='kegg diamond比对阈值')
    parser.add_argument('--kegg_ko_gene_relation', help='自建库, ko与acc')
    parser.add_argument('--kegg_relation_databse', help='自建数据库, kegg各个item间关系')
    parser.add_argument('--cazy_diamond_database', help='cazy 蛋白序列数据库')
    parser.add_argument('--cazy_diamond_evalue', help='cazy diamond比对阈值')
    parser.add_argument('--cazy_relation_database', help='cazy 各个item间关系')
    parser.add_argument('--eggnog_diamond_database', help='eggnog 蛋白序列数据库')
    parser.add_argument('--eggnog_diamond_evalue', help='eggnog diamond 比对阈值')
    parser.add_argument('--eggnog_relation_database', help='eggnog 各个item间关系')
    ## 报告参数 
    parser.add_argument('--report', action='store_true', default=False)
    parser.add_argument('--report_pdf', help='是否生成pdf报告')
    parser.add_argument('--report_english', help='是否生成英文版报告', default=False) 

    args = vars(parser.parse_args())
    MainPipeline(args).start()
