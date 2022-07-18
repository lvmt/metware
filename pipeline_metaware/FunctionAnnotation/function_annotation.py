#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-30 09:05:22
'''


'''
功能注释
'''

from ctypes import util
import sys 
import os 
import textwrap

try:
    from Utils import utils
except Exception as e:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
    from Utils import utils


class Function:

    def __init__(self, args):
        self.args = args 


    def start(self):
        KEGG(self.args).start()
        CAZy(self.args).start()
        eggNOG(self.args).start()



class KEGG:
    # KEGG 功能注释 
    def __init__(self, args):
        self.args = args
        self.projdir = self.args['projdir']
        self.sample_file = self.args['sample_file']
        self.kegg_diamond_evalue = self.args['kegg_diamond_evalue']
        self.kegg_diamond_database = self.args['kegg_diamond_database']
        self.kegg_ko_gene_relation = self.args['kegg_ko_gene_relation']
        self.kegg_relation_database = self.args['kegg_relation_database']
        self.kegg_databse = self.args['kegg_database']
        self.pathway_description = self.args['pathway_description']


    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/KEGG')
        cmd = textwrap.dedent(f'''
        echo "start time: "  `date`
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d {self.kegg_diamond_database} \\
            -q {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.kegg_diamond_evalue} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore &&\\
        
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/blastp.sh'
        utils.write_cmd(cmd, shellname)

    
    def anno(self):
        # kegg数据库注释 
        cmd = textwrap.dedent(f'''
        echo "start time: "  `date`
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8 \\
            --ko_gene {self.kegg_ko_gene_relation} \\
            --ko_relation {self.kegg_relation_database} \\
            --result {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8.anno &&\\

        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/anno.sh'
        utils.write_cmd(cmd, shellname)
        
    
    def abundance_stat(self):
        # 给kegg注释结果增加物种和丰度信息 
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        ## kegg 注释结果给予预测基因,添加物种注释信息和丰度信息, 并拆分至各个层级
        mkdir -p {self.projdir}/5.FunctionAnnotation/KEGG/Relative
        mkdir -p {self.projdir}/5.FunctionAnnotation/KEGG/Absolute

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_stat_abundance.py \\
            --kegg_anno {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8.anno \\
            --abundance_table {self.projdir}/4.TaxAnnotation/Absolute/Unigenes.tax_abundance.absolute.total.xls \\
            --sample_file {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.kegg \\
            --result_dir {self.projdir}/5.FunctionAnnotation/KEGG &&\\

        ## 绘制pathway通路图
        ~/.conda/envs/python3_lmt/bin/python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/plot_kegg_pathway.py \\
            --kegg_tax {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8.anno.kegg.tax.xls \\
            --kegg_database {self.kegg_databse} \\
            --pathway_description {self.pathway_description} \\
            --result_dir {self.projdir}/5.FunctionAnnotation/KEGG/pathway_map &&\\

        ## 生成pathway 报告 
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/pathway_map_report.py \\
            --projdir {self.projdir} \\
            --stat_info {self.projdir}/5.FunctionAnnotation/KEGG/pathway_map/stat.xls \\
            --template_dir ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/templates \\
            --report_dir {self.projdir}/5.FunctionAnnotation/KEGG/pathway_report/  &&\\

        echo "end time: "  `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)


    def gene_state(self):
        # 统计不同kegg levels的基因计数
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        ## 不同水平的基因数目统计
        ## 汇总统计|单样本统计
        mkdir -p {self.projdir}/5.FunctionAnnotation/KEGG/GeneStat
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_stat_gene.py \\
            --kegg_anno {self.projdir}/5.FunctionAnnotation/KEGG/Unigenes.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.gene.readsNum \\
            --sample_file {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.kegg.gene_stat  \\
            --result_dir {self.projdir}/5.FunctionAnnotation/KEGG &&\\

        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/gene_stat.sh'
        utils.write_cmd(cmd, shellname)
         

    def start(self):
        self.blastp()
        self.anno()
        self.abundance_stat()
        self.gene_state()



class CAZy:
    def __init__(self, args):
        self.args = args
        self.projdir = self.args['projdir']
        self.sample_file = self.args['sample_file']
        self.cazy_diamond_databse = self.args['cazy_diamond_databse']
        self.cazy_diamond_evalue = self.args['cazy_diamond_evalue']
        self.cazy_relation_database = self.args['cazy_relation_database']


    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/CAZy')
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d {self.cazy_diamond_databse} \\
            -q {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/CAZy/Unigenes.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.cazy_diamond_evalue} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore &&\\
        
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/blastp.sh'
        utils.write_cmd(cmd, shellname)
    

    def anno(self):
        # kegg数据库注释 
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/CAZy/Unigenes.protein.m8 \\
            --cazy_relation {self.cazy_relation_database} \\
            --result {self.projdir}/5.FunctionAnnotation/CAZy/Unigenes.protein.m8.anno &&\\

        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/anno.sh'
        utils.write_cmd(cmd, shellname)

    
    def abundance_stat(self):
        # 给cazy注释结果增加物种和丰度信息 
        cmd = textwrap.dedent(f'''
        echo "start time: "  `date`
        mkdir -p {self.projdir}/5.FunctionAnnotation/CAZy/Relative
        mkdir -p {self.projdir}/5.FunctionAnnotation/CAZy/Absolute

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_stat_abundance.py \\
            --cazy_anno {self.projdir}/5.FunctionAnnotation/CAZy/Unigenes.protein.m8.anno  \\
            --abundance_table {self.projdir}/4.TaxAnnotation/Absolute/Unigenes.tax_abundance.absolute.total.xls \\
            --sample_file  {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.cazy \\
            --result_dir {self.projdir}/5.FunctionAnnotation/CAZy &&\\
        
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)


    def gene_stat(self):
        # 统计不同cazy水平的基因计数
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        mkdir -p {self.projdir}/5.FunctionAnnotation/CAZy/GeneStat
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_stat_gene.py \\
            --cazy_anno {self.projdir}/5.FunctionAnnotation/CAZy/Unigenes.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.gene.readsNum  \\
            --sample_file {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.cazy.gene_stat  \\
            --result_dir {self.projdir}/5.FunctionAnnotation/CAZy &&\\

        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/gene_stat.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        self.blastp()
        self.anno()
        self.abundance_stat()
        self.gene_stat()



class eggNOG:
    
    def __init__(self, args):
        self.args = args
        self.projdir = self.args['projdir']
        self.sample_file = self.args['sample_file']
        self.eggnog_diamond_evalue = self.args['eggnog_diamond_evalue']
        self.eggnog_diamond_databse = self.args['eggnog_diamond_databse']
        self.eggnog_relation_database = self.args['eggnog_relation_database']


    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/eggNOG')
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d  {self.eggnog_diamond_databse} \\
            -q {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/eggNOG/Unigenes.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.eggnog_diamond_evalue} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore &&\\
        
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/blastp.sh'
        utils.write_cmd(cmd, shellname)


    def anno(self):
        cmd = textwrap.dedent(f'''
        echo  "start time: " `date`
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/eggNOG/Unigenes.protein.m8 \\
            --eggNOG_relation {self.eggnog_relation_database} \\
            --result {self.projdir}/5.FunctionAnnotation/eggNOG/Unigenes.protein.m8.anno &&\\
        
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/anno.sh'
        utils.write_cmd(cmd, shellname)


    def abundance_stat(self):
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        mkdir -p {self.projdir}/5.FunctionAnnotation/eggNOG/Relative
        mkdir -p {self.projdir}/5.FunctionAnnotation/eggNOG/Absolute

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_stat_abundance.py \\
            --eggNOG_anno {self.projdir}/5.FunctionAnnotation/eggNOG/Unigenes.protein.m8.anno \\
            --abundance_table {self.projdir}/4.TaxAnnotation/Absolute/Unigenes.tax_abundance.absolute.total.xls \\
            --sample_file {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.eggNOG \\
            --result_dir {self.projdir}/5.FunctionAnnotation/eggNOG/ &&\\

        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)
         

    def gene_stat(self):
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        mkdir -p {self.projdir}/5.FunctionAnnotation/eggNOG/GeneStat 
        
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_stat_gene.py  \\
            --eggNOG_anno {self.projdir}/5.FunctionAnnotation/eggNOG/Unigenes.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/Unigenes.total.gene.readsNum \\
            --sample_file {self.sample_file} \\
            --result_suffix Unigenes.protein.m8.anno.eggNOG.gene_stat \\
            --result_dir {self.projdir}/5.FunctionAnnotation/eggNOG &&\\
            
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/gene_stat.sh'
        utils.write_cmd(cmd, shellname)


    def start(self):
        self.blastp()
        self.anno()
        self.abundance_stat()
        self.gene_stat()

        


if __name__ == '__main__':
    print('单例测试')
    

