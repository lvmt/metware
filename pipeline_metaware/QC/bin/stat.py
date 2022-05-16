#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-16 11:11:13
'''

# 质控图片绘制
import sys  
import json 
import pandas as pd 
import matplotlib.pyplot as plt
plt.switch_backend('agg')


json_file = sys.argv[1]
out_suffix = json_file.replace('.fastp.json', '')

info = json.load(open(json_file, 'r'))
read1_quality_df = pd.DataFrame.from_dict(info['read1_after_filtering']['quality_curves'])
read1_GC_df = pd.DataFrame.from_dict(info['read1_after_filtering']['content_curves'])
read2_quality_df = pd.DataFrame.from_dict(info['read2_after_filtering']['quality_curves'])
read2_GC_df = pd.DataFrame.from_dict(info['read2_after_filtering']['content_curves'])


# 碱基质量值绘图
fig1_qua = read1_quality_df['mean'].plot()
plt.ylim(0,40)
plt.xlabel('position')
plt.ylabel('quality')
plt.title('read1')
fig1_qua.figure.savefig(f'{out_suffix}.read1.qual.png', dpi=600)

fig2_qua = read2_quality_df['mean'].plot()
plt.ylim(0,40)
plt.xlabel('position')
plt.ylabel('quality')
plt.title('read2')
fig2_qua.figure.savefig(f'{out_suffix}.read2.qual.png', dpi=600)


# GC含量绘图
fig1_gc = read1_GC_df.plot()
plt.ylim(0.1, 0.6)
plt.xlabel('position')
plt.ylabel('percent')
plt.title('read1')
fig1_gc.figure.savefig(f'{out_suffix}.read1.GC.png', dpi=600)

fig2_gc = read1_GC_df.plot()
plt.ylim(0.1, 0.6)
plt.xlabel('position')
plt.ylabel('percent')
plt.title('read2')
fig1_gc.figure.savefig(f'{out_suffix}.read2.GC.png', dpi=600)

# 生成统计文件
before_df = pd.DataFrame([info['summary']['before_filtering']]).rename(index={0:'before'})
after_df = pd.DataFrame([info['summary']['after_filtering']]).rename(index={0:'after'})
combine_df = before_df.append(after_df)
combine_df.to_csv(f'{out_suffix}.stat.xls', sep='\t')




