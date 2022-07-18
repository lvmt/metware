import os
import re
import sys
import pypandoc
import argparse
import pandas as pd
import glob 
from jinja2 import Environment, FileSystemLoader



REPORT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(REPORT_DIR))
sys.path.append(ROOT_DIR)


## paras define
parser = argparse.ArgumentParser(description='报告解读')
parser.add_argument('--projdir', help='项目分析目录')
parser.add_argument('--sample_file', help='样本信息配置文件')
parser.add_argument('--template_dir', help='报告模板')
parser.add_argument('--report', help='结题报告输出名称')
parser.add_argument('--pdf', help='是否输出pdf结果', choices=['yes', 'no'], default='no')
parser.add_argument('--project_num', help='项目编号 XXXXXX-12138', default='XXXXXX-12138')
parser.add_argument('--project_name', help='项目名称, xxxxxxx-250941', default='xxxxxxx-250941')
parser.add_argument('--report_time', help='报告时间 2022-02-20', default='2022-02-20')
parser.add_argument('--report_num', help='报告编号: xxxxxx-54188', default='xxxxxx-54188')


args = vars(parser.parse_args())


def get_samples(sample_file):
    samples = []
    with open(sample_file, 'r') as fr:
        for line in fr:
            if line.startswith('#'):
                continue
            linelist = line.strip('\n').split('\t')
            samples.append(linelist[0])
    return samples 

samples = get_samples(args['sample_file'])
############################################################################

## 获取模板变量
templates_dir = args.get('template_dir')
env = Environment(loader=FileSystemLoader(templates_dir))
template = env.get_template('report_met.html')
docx_template = env.get_template('report_met.v2_docx.html')


##################################### 解析文档填充内容 
## 样本基本信息
sample_info = []
with open(args['sample_file'], 'r') as fr:
    for line in fr:
        if line.startswith('#'):
            continue
        linelist = line.strip('\n').split('\t')
        sampleid = linelist[0]
        group = linelist[1]
        sample_info.append([sampleid, group])


###############################  qc #################################################
os.system('mkdir -p {report}/src/pictures/1.Clean/'.format(**args))
qc_stat_list = []
qc_stat_file = '{projdir}/1.Clean/all.samples.stat.txt'.format(**args)
with open(qc_stat_file, 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        qc_stat_list.append(linelist)

gc_png_list = []
for sample in samples:
    os.system('cp {projdir}/1.Clean/{sample}/{sample}.read1.qual.png {report}/src/pictures/1.Clean/'.format(**locals(), **args))
    os.system('cp {projdir}/1.Clean/{sample}/{sample}.read2.qual.png {report}/src/pictures/1.Clean/'.format(**locals(), **args))
    os.system('cp {projdir}/1.Clean/{sample}/{sample}.read1.GC.png {report}/src/pictures/1.Clean/'.format(**locals(), **args))
    os.system('cp {projdir}/1.Clean/{sample}/{sample}.read2.GC.png {report}/src/pictures/1.Clean/'.format(**locals(), **args))
    gc_png_list.extend((f'src/pictures/1.Clean/{sample}.read1.GC.png', f'src/pictures/1.Clean/{sample}.read2.GC.png'))


########################################### assembly ###################################################
os.system('mkdir -p {report}/src/pictures/2.Assembly/'.format(**args))
os.system('cp {projdir}/2.Assembly/all.samples.assembly_assessment.txt {report}/src/pictures/2.Assembly/'.format(**args))
assembly_stat_list = []
assembly_stat_file = '{report}/src/pictures/2.Assembly/all.samples.assembly_assessment.txt'.format(**args)
ass_df = pd.read_csv(assembly_stat_file, sep='\t', index_col=0).T.reset_index()
assembly_stat_list.append(ass_df.columns)
assembly_stat_list.extend(map(list, ass_df.values))


########################################## gene prediction #############################################
os.system('mkdir -p {report}/src/pictures/3.predict/'.format(**args))
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.total.nucleotide.summary.txt {report}/src/pictures/3.predict/'.format(**args))
uniq_gene_stat_list = []
uniq_gene_stat_file = '{report}/src/pictures/3.predict/Unigenes.total.nucleotide.summary.txt'.format(**args)
with open(uniq_gene_stat_file, 'r') as fr:
    for line in fr:
        linelist = line.strip('\n').split('\t')
        uniq_gene_stat_list.append(linelist)

os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.total.nucleotide.fa.len.png {report}/src/pictures/3.predict/'.format(**args))  
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.plot.core.png {report}/src/pictures/3.predict/'.format(**args))  
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.plot.pan.png {report}/src/pictures/3.predict/'.format(**args))   
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.plot.boxplot.png {report}/src/pictures/3.predict/'.format(**args))   
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.plot.veen.png {report}/src/pictures/3.predict/'.format(**args))   
os.system('cp {projdir}/3.GenePrediction/Cluster/Unigenes.plot.cor.png {report}/src/pictures/3.predict/'.format(**args))   


########################################## Taxonomy annotation #############################################
os.system('mkdir -p {report}/src/pictures/4.Taxonomy/'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/Krona/krona.html {report}/src/pictures/4.Taxonomy/'.format(**args))
#  门水平和属水平的物种相对丰度柱形图（样品）
os.system('cp {projdir}/4.TaxAnnotation/Plot/top/top10.phylum.barplot.png  {report}/src/pictures/4.Taxonomy/p_top10_sample.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/top/top10.genus.barplot.png  {report}/src/pictures/4.Taxonomy/g_top10_sample.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/top_group/top10.phylum.barplot.png  {report}/src/pictures/4.Taxonomy/p_top10_group.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/top_group/top10.genus.barplot.png  {report}/src/pictures/4.Taxonomy/g_top10_group.png'.format(**args))
# 热图
os.system('cp {projdir}/4.TaxAnnotation/Plot/GeneNums.BetweenSamples.heatmap/taxonomy.phylum.heatmap.png  {report}/src/pictures/4.Taxonomy/gene_nums_heatmap.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/heatmap/taxonomy.phylum.heatmap.png  {report}/src/pictures/4.Taxonomy/relative_abundance_heatmap.png'.format(**args))
# 基于物种丰度降维分析
os.system('cp {projdir}/4.TaxAnnotation/Plot/PCA/phylum_PCA.png  {report}/src/pictures/4.Taxonomy/p_pca.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/NMDS/phylum_NMDS.png  {report}/src/pictures/4.Taxonomy/p_nmds.png'.format(**args))
os.system('cp {projdir}/4.TaxAnnotation/Plot/PCoA/phylum_PCoA.cluster.png  {report}/src/pictures/4.Taxonomy/p_pcoa.png'.format(**args))
# anosim
anosim_png = glob.glob('{projdir}/4.TaxAnnotation/Plot/Anosim/phylum/*png'.format(**args))[0]
os.system('cp {anosim_png}  {report}/src/pictures/4.Taxonomy/p_anosim.png'.format(**args, **locals()))
# 基于物种丰度的样品聚类分析
os.system('cp {projdir}/4.TaxAnnotation/Plot/UPGMA/bar_tree_p10.png  {report}/src/pictures/4.Taxonomy/bar_tree_p10.png'.format(**args))
# metastat 
metastat_png = glob.glob('{projdir}/4.TaxAnnotation/Plot/Metastats/*/*.png'.format(**args))[0]
os.system('cp {metastat_png}  {report}/src/pictures/4.Taxonomy/metastat.png'.format(**args, **locals()))
# LEfse
lda_png = glob.glob('{projdir}/4.TaxAnnotation/Plot/LEfSe/*/*LDA.png'.format(**args))[0]
lda_tree_png = glob.glob('{projdir}/4.TaxAnnotation/Plot/LEfSe/*/*LDA.tree.png'.format(**args))[0]
os.system('cp {lda_png}  {report}/src/pictures/4.Taxonomy/LDA.png'.format(**args, **locals()))
os.system('cp {lda_tree_png}  {report}/src/pictures/4.Taxonomy/LDA.tree.png'.format(**args, **locals()))


########################################## Taxonomy function #############################################
os.system('mkdir -p {report}/src/pictures/5.FunctionAnnotation/'.format(**args))
# 注释基因数目统计
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.unigenes.num.png {report}/src/pictures/5.FunctionAnnotation/kegg.unigenes.num.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.unigenes.num.png {report}/src/pictures/5.FunctionAnnotation/eggnog.unigenes.num.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.unigenes.num.png {report}/src/pictures/5.FunctionAnnotation/cazy.unigenes.num.png'.format(**args))
func_gene_count_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_gene_count_png_list.append('src/pictures/5.FunctionAnnotation/{_class}.unigenes.num.png'.format(**args, **locals()))

# 功能相对丰度概况
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.unigenes.level1.bar.png {report}/src/pictures/5.FunctionAnnotation/kegg.unigenes.level1.bar.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.unigenes.level1.bar.png {report}/src/pictures/5.FunctionAnnotation/eggnog.unigenes.level1.bar.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.unigenes.level1.bar.png {report}/src/pictures/5.FunctionAnnotation/cazy.unigenes.level1.bar.png'.format(**args))
func_relative_abundance_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_relative_abundance_png_list.append('src/pictures/5.FunctionAnnotation/{_class}.unigenes.level1.bar.png'.format(**locals(), **args))

# 功能相对丰度聚类分析
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.cluster.ko.png {report}/src/pictures/5.FunctionAnnotation/kegg.cluster.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.cluster.og.png {report}/src/pictures/5.FunctionAnnotation/eggnog.cluster.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.cluster.level2.png {report}/src/pictures/5.FunctionAnnotation/cazy.cluster.png'.format(**args))
func_abundance_cluster_heatmap_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_abundance_cluster_heatmap_png_list.append('src/pictures/5.FunctionAnnotation/{_class}.cluster.png'.format(**locals(), **args))

# 基于功能丰度的降维分析
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.pca.png {report}/src/pictures/5.FunctionAnnotation/kegg.pca.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.pca.png {report}/src/pictures/5.FunctionAnnotation/eggnog.pca.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.pca.png {report}/src/pictures/5.FunctionAnnotation/cazy.pca.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.nmds.png {report}/src/pictures/5.FunctionAnnotation/kegg.nmds.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.nmds.png {report}/src/pictures/5.FunctionAnnotation/eggnog.nmds.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.nmds.png {report}/src/pictures/5.FunctionAnnotation/cazy.nmds.png'.format(**args))
func_abundance_pca_nmds_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_abundance_pca_nmds_png_list.append(['src/pictures/5.FunctionAnnotation/{_class}.pca.png'.format(**args, **locals()), 'src/pictures/5.FunctionAnnotation/{_class}.nmds.png'.format(**locals(), **args)])

# 于功能丰度的Bray-Curtis 距离的降维分析
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.pcoa.png {report}/src/pictures/5.FunctionAnnotation/kegg.pcoa.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.pcoa.png {report}/src/pictures/5.FunctionAnnotation/eggnog.pcoa.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.pcoa.png {report}/src/pictures/5.FunctionAnnotation/cazy.pcoa.png'.format(**args))
func_abundance_pcoa_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_abundance_pcoa_png_list.append(f'src/pictures/5.FunctionAnnotation/{_class}.pcoa.png')

# 基于功能丰度的Anosim分析
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.anosim.png {report}/src/pictures/5.FunctionAnnotation/kegg.anosim.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.anosim.png {report}/src/pictures/5.FunctionAnnotation/eggnog.anosim.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.anosim.png {report}/src/pictures/5.FunctionAnnotation/cazy.anosim.png'.format(**args))
func_abundance_anosim_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_abundance_anosim_png_list.append('src/pictures/5.FunctionAnnotation/{_class}.anosim.png'.format(**locals(), **args))

# 基于功能丰度的样品聚类分析
os.system('cp {projdir}/5.FunctionAnnotation/KEGG/Plot/kegg.unigenes.level1.bar.tree.png {report}/src/pictures/5.FunctionAnnotation/kegg.unigenes.level1.bar.tree.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/eggNOG/Plot/eggnog.unigenes.level1.bar.tree.png {report}/src/pictures/5.FunctionAnnotation/eggnog.unigenes.level1.bar.tree.png'.format(**args))
os.system('cp {projdir}/5.FunctionAnnotation/CAZy/Plot/cazy.unigenes.level1.bar.tree.png {report}/src/pictures/5.FunctionAnnotation/cazy.unigenes.level1.bar.tree.png'.format(**args))
func_abundance_cluster_tree_png_list = []
for _class in ['kegg', 'eggnog', 'cazy']:
    func_abundance_cluster_tree_png_list.append('src/pictures/5.FunctionAnnotation/{_class}.unigenes.level1.bar.tree.png'.format(**locals(), **args))

# pathway_map
os.system('cp -r {projdir}/5.FunctionAnnotation/KEGG/pathway_report {report}/src/pictures/5.FunctionAnnotation/pathway_report'.format(**args))


# 组间功能差异的Metastat分析, 可能有结果，可能没有结果
func_metastat_png_list = []
metastat_png = glob.glob('{projdir}/5.FunctionAnnotation/*/Plot/*.metastat.png'.format(**args))
for metastat_png_single in metastat_png:
    os.system('cp {metastat_png_single} {report}/src/pictures/5.FunctionAnnotation/'.format(**locals(), **args))
    metastat_png_single = os.path.basename(metastat_png_single)
    func_metastat_png_list.append(f'src/pictures/5.FunctionAnnotation/{metastat_png_single}')

# 根据组间具有差异的功能进行PCA分析并绘制丰度聚类热图
func_metastat_pca_heatmap_png_list = []
diff_fun_pca_pngs = glob.glob('{projdir}/5.FunctionAnnotation/*/Plot/*.diff_func.pca.png'.format(**args))
diff_fun_heatmap_pngs = glob.glob('{projdir}/5.FunctionAnnotation/*/Plot/*.diff_func.heatmap.png'.format(**args))
for pca,heatmap in zip(diff_fun_pca_pngs, diff_fun_heatmap_pngs):
    os.system('cp {pca} {report}/src/pictures/5.FunctionAnnotation/'.format(**args, **locals()))
    os.system('cp {heatmap} {report}/src/pictures/5.FunctionAnnotation/'.format(**args, **locals()))
    pca = os.path.basename(pca)
    heatmap = os.path.basename(heatmap)
    func_metastat_pca_heatmap_png_list.append([f'src/pictures/5.FunctionAnnotation/{pca}', f'src/pictures/5.FunctionAnnotation/{heatmap}'])

############################################### 抗性基因注释 #############################################################
os.system('mkdir -p {report}/src/pictures/6.CARD/'.format(**args))
# 抗性基因丰度概况
os.system('cp {projdir}/6.CARD/Plot/card.relative_abundance_top20.png {report}/src/pictures/6.CARD/card.relative_abundance_top20.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/card.relative_abundance_ratio_top20.png {report}/src/pictures/6.CARD/card.relative_abundance_ratio_top20.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/aro.top10.circos.png {report}/src/pictures/6.CARD/aro.top10.circos.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/card.mechanism.taxonomy.png {report}/src/pictures/6.CARD/card.mechanism.taxonomy.png'.format(**args))
# 抗性基因类型的分布及其丰度聚类分析
os.system('cp {projdir}/6.CARD/Plot/aro_heatmap_black.png {report}/src/pictures/6.CARD/aro_heatmap_black.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/aro_heatmap_white.png {report}/src/pictures/6.CARD/aro_heatmap_white.png'.format(**args))
# 基于抗性基因丰度的Anosim分析
os.system('cp {projdir}/6.CARD/Plot/card_anosim_show.png {report}/src/pictures/6.CARD/card_anosim_show.png'.format(**args))
#  组间抗性基因数目差异分析
os.system('cp {projdir}/6.CARD/Plot/gene_nums_box.png {report}/src/pictures/6.CARD/gene_nums_box.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/aro_nums_box.png {report}/src/pictures/6.CARD/aro_nums_box.png'.format(**args))
os.system('cp {projdir}/6.CARD/Plot/card_venn_flower_display.png {report}/src/pictures/6.CARD/card_venn_flower_display.png'.format(**args))



############################# 模板填充 ###############################################
os.system('mkdir -p {report}/src'.format(**args))
os.system('cp -r templates/src/* {report}/src'.format(**args))
outhtml = os.path.abspath(os.path.join(args['report'], 'report.final.html'))


context = {}
context['project_num'] = args['project_num']
context['project_name'] = args['project_name']
context['report_num'] = args['report_num']
context['report_time'] = args['report_time']
context['sample_info'] = sample_info
context['qc_stat_list'] = qc_stat_list
context['gc_png_list'] = gc_png_list
context['assembly_stat_list'] = assembly_stat_list
context['uniq_gene_stat_list'] = uniq_gene_stat_list

# 更新功能注释参数
context['func_gene_count_png_list'] = func_gene_count_png_list
context['func_relative_abundance_png_list'] = func_relative_abundance_png_list
context['func_abundance_cluster_heatmap_png_list'] = func_abundance_cluster_heatmap_png_list
context['func_abundance_pca_nmds_png_list'] = func_abundance_pca_nmds_png_list
context['func_abundance_pcoa_png_list'] = func_abundance_pcoa_png_list
context['func_abundance_anosim_png_list'] = func_abundance_anosim_png_list
context['func_abundance_cluster_tree_png_list'] = func_abundance_cluster_tree_png_list
context['func_metastat_png_list'] = func_metastat_png_list
context['func_metastat_pca_heatmap_png_list'] = func_metastat_pca_heatmap_png_list



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
        ~/pipeline/metagenomics/software/wkhtmltox/bin/wkhtmltopdf \\
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
# 由于在html中, 部分表格存在筛选框, 需要利用正则表达式将html中相关的内容去掉
docx_content = docx_template.render(context)
docx_html = outhtml.replace('html', 'docx.html')
with open(docx_html, 'w') as fw:
    fw.write(docx_content)

os.chdir(os.path.abspath(os.path.dirname(docx_html)))
docx_result = outhtml.replace('html', 'docx')
pypandoc.convert_file(docx_html, 'docx', outputfile=docx_result)
os.system(f'rm -f {docx_html}')



 