import argparse
import pandas as pd
import pydoc
from jinja2 import Environment, FileSystemLoader
import os
import re
import sys
import pypandoc



REPORT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(REPORT_DIR))
sys.path.append(ROOT_DIR)


## paras define
parser = argparse.ArgumentParser(description='报告解读')
parser.add_argument('--projdir', help='项目分析目录')
parser.add_argument('--sample_file', help='样本信息配置文件')
parser.add_argument('--template_dir', help='报告模板')
parser.add_argument('--report', help='截图报告输出')
parser.add_argument('--pdf', help='是否输出pdf结果', default='yes')


args = vars(parser.parse_args())


## 获取模板变量
templates_dir = args.get('template_dir')
env = Environment(loader=FileSystemLoader(templates_dir))
template = env.get_template('report_met.html')


######################## 解析文档填充内容 

## qc 
qc_stat_list = []
qc_stat_file = '{projdir}/1.Clean/all.samples.stat.txt'.format(**args)
with open(qc_stat_file, 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        qc_stat_list.append(linelist)


## assembly
assembly_stat_list = []
assembly_stat_file = '{projdir}/2.Assembly/all.samples.assembly_assessment.txt'.format(**args)
ass_df = pd.read_csv(assembly_stat_file, sep='\t', index_col=0).T.reset_index()
assembly_stat_list.append(ass_df.columns)
assembly_stat_list.extend(map(list, ass_df.values))
print(assembly_stat_list)
# with open(assembly_stat_file, 'r') as fr:
#     for line in fr:
#         linelist = line.strip('\n').split('\t')
#         assembly_stat_list.append(linelist)


## gene prediction 
uniq_gene_stat_list = []
uniq_gene_stat_file = '{projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.summary.txt'.format(**args)
with open(uniq_gene_stat_file, 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        uniq_gene_stat_list.append(linelist)
os.system('mkdir -p {report}/src/pictures/3.predict/'.format(**args))
os.system('cp {projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa.len.png {report}/src/pictures/3.predict/'.format(**args))





############################# 模板填充 
os.system('mkdir -p {report}/src'.format(**args))
os.system('cp -r templates/src/* {report}/src'.format(**args))
outhtml = os.path.abspath(os.path.join(args['report'], 'report.final.html'))


context = {}
context['qc_stat_list'] = qc_stat_list
context['assembly_stat_list'] = assembly_stat_list
context['uniq_gene_stat_list'] = uniq_gene_stat_list


## 输出pdf报告
pdf_report = args['pdf']
if pdf_report == 'yes':
    temp_html = outhtml.replace('html', 'temp.html')
    with open(temp_html, 'w') as temp:
        context['pdf_report'] = True
        pdf_html = template.render(context)
        temp.write(pdf_html)

    pdf_result = outhtml.replace('html', 'pdf')
    cmd = '''
        wkhtmltopdf \\
            -n \\
            --print-media-type \\
            --images \\
            --page-size A4 \\
            --debug-javascript \\
            --enable-javascript \\
            --no-stop-slow-scripts \\
            --javascript-delay 15000 \\
            {temp_html} \\
            {pdf_result} 
            
        '''.format(**locals())
    os.system(cmd)
    os.system(f'rm -f {temp_html}')


## 输出html报告
content = template.render(context)
with open(outhtml, 'w') as fw:
    fw.write(content)


## 输出docx报告
os.chdir(os.path.abspath(os.path.dirname(outhtml)))
docx_result = outhtml.replace('html', 'docx')
pypandoc.convert_file(outhtml, 'docx', outputfile=docx_result)








 