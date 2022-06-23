#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-17 16:15:19
'''

# 生成基因预测结果的汇总信息
import os 
import pandas as pd  
import matplotlib.pyplot as plt
plt.switch_backend('agg')



def main(args):
    result_suffix = args['result_suffix']
    name = os.path.basename(args['infile']).split('.')[0]
    df = pd.read_csv(args['infile'], sep='\t')
    orf_nums = df.shape[0]
    total_len = round(df['Length'].sum() / 1000 / 1000, 2)  # M
    both_yes = df[(df['Initiation_codon'] == 'yes') & (df['Stop_codon']=='yes')].shape[0]
    just_start = df[(df['Initiation_codon'] == 'yes') & (df['Stop_codon']=='no')].shape[0]
    just_end = df[(df['Initiation_codon'] == 'no') & (df['Stop_codon']=='yes')].shape[0]
    both_no = df[(df['Initiation_codon'] == 'no') & (df['Stop_codon']=='no')].shape[0]
    average_len = round(df['Length'].mean(), 2)
    gc_percent = round(df['GC(%)'].mean(), 2)
    both_yes_per = round(both_yes / orf_nums * 100, 2)
    just_start_per = round(just_start / orf_nums * 100, 2)
    just_end_per = round(just_end / orf_nums * 100, 2)
    both_no_per = round(both_no / orf_nums * 100, 2)

    info = {
        'sample': name,
        'ORFs NO.': orf_nums,
        'integrity:all': f'{both_yes}({both_yes_per}%)',
        'integrity:none': f'{both_no}({both_no_per}%)',
        'integrity:end': f'{just_end}({just_end_per}%)',
        'integrity:start': f'{just_start}({just_start_per}%)',
        'Total Len.(Mbp)': total_len,
        'Average Len.(bp)': average_len,
        'GC percent': gc_percent
    }   
    
    outfile = f'{result_suffix}.summary.txt'
    with open(outfile, 'w') as fw:
        for key,value in info.items():
            fw.write(f'{key}\t{value}\n')

    ## 绘制长度分箱计数分布图
    bins = list(range(100, 1700, 100))  # 长度从100开始,因为前面已经过滤掉了<100nt的预测结果
    labels = list(map(str, bins[:-1]))
    df['bin'] = pd.cut(df.Length, bins=bins, labels=labels, include_lowest=True, right=False) # 左闭右开
    df['bin'] = df['bin'].fillna(labels[-1])  # 分箱时,对于>1600nt,被赋值为NaN; 此举表示最后一个分箱为>=1500nt
    x =  list(df['bin'].value_counts().sort_index().index)
    y =  list(df['bin'].value_counts().sort_index().values)
    plt.plot(x, y)
    plt.xticks(rotation=270)
    plt.xlabel('ORF Length(nt)')
    plt.ylabel('Counts')
    plt.title(f'{name} Predict Gene Length Distribution')
    plt.savefig(f'{result_suffix}.fa.len.png', dpi=400, bbox_inches='tight')
    


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='基于统计信息生成汇总信息')
    parser.add_argument('--infile', help='输入文件: 基因信息统计文件')
    parser.add_argument('--result_suffix', help='输出前缀')


    args = vars(parser.parse_args())
    main(args)



