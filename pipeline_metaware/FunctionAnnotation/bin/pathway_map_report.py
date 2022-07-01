#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-30 16:43:36
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-30 16:43:36
'''


import os  
import sys 
import argparse
from jinja2 import Environment, FileSystemLoader

REPORT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(REPORT_DIR))
sys.path.append(ROOT_DIR)


## paras define
parser = argparse.ArgumentParser(description='报告解读')
parser.add_argument('--projdir', help='项目分析目录')
parser.add_argument('--stat_info', help='统计信息')
parser.add_argument('--template_dir', help='报告模板')
parser.add_argument('--report_dir', help='报告输出')

args = vars(parser.parse_args())


## 获取模板变量
templates_dir = args.get('template_dir')
env = Environment(loader=FileSystemLoader(templates_dir))
template = env.get_template('pathway_map_report.html') 


def get_pathway_map_list(stat_info):
    pathway_map_list = []
    with open(stat_info, 'r') as fr:
        for line in fr:
            if not line.startswith(('map', 'ko')):
                continue
            linelist = line.strip('\n').split('\t')
            _map = linelist[0]
            linelist[0] = '<a href="src/{_map}.png" target="_blank">{_map}</a>'.format(**locals())
            pathway_map_list.append(linelist)

    return pathway_map_list 


pathway_map_list = get_pathway_map_list(args['stat_info'])

####################################################
context = {}
context['pathway_map_list'] = pathway_map_list


## 输出html报告
os.system('mkdir -p {report_dir}/src'.format(**args))
os.system('ln -sf {projdir}/5.FunctionAnnotation/KEGG/pathway_map/map*png {report_dir}/src'.format(**args))
outhtml = '{report_dir}/pathway_map.html'.format(**args)
content = template.render(context)
with open(outhtml, 'w') as fw:
    fw.write(content)

