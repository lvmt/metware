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



class KEGG:
    # KEGG 功能注释 
    def __init__(self, args):
        self.args = args
        self.db = self.args['db_kegg']
        self.projdir = self.args['projdir']
        self.threshold = self.args['threshold']


    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/KEGG')
        cmd = textwrap.dedent(f'''
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d ~/pipeline/metagenomics/database/kegg/KEGG.protein.seq.microbe.dmnd \\
            -q {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.threshold} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/blastp.sh'
        utils.write_cmd(cmd, shellname)

    
    def anno(self):
        # kegg数据库注释 
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8 \\
            --ko_gene ~/pipeline/metagenomics/database/kegg/KO_gene \\
            --ko_relation ~/pipeline/metagenomics/database/kegg/kegg_relation.json \\
            --result {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8.anno 
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/anno.sh'
        utils.write_cmd(cmd, shellname)
        
    
    def abundance_stat(self):
        # 给kegg注释结果增加物种和丰度信息 
        cmd = textwrap.dedent(f'''
        ## kegg 注释结果给予预测基因,添加物种注释信息和丰度信息, 并拆分至各个层级

        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_stat_abundance.py \\
            --kegg_anno {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8.anno \\
            --abundance_table {self.projdir}/4.TaxAnnotation/NonRundant.tax_abundance.absolute.total.xls \\
            --sample_file ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt \\
            --result_suffix {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8.anno.kegg

        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/KEGG/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)


    def gene_state(self):
        # 统计不同kegg levels的基因计数
        cmd = textwrap.dedent(f'''
        ## 不同水平的基因数目统计
        ## 汇总统计|单样本统计
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/kegg_stat_gene.py \\
            --kegg_anno {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum \\
            --sample_file ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt \\
            --result_suffix {self.projdir}/5.FunctionAnnotation/KEGG/NonRundant.protein.m8.anno.kegg.gene_stat 
        
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
        self.threshold = self.args['threshold']


    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/CAZy')
        cmd = textwrap.dedent(f'''
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d ~/pipeline/metagenomics/database/CAZy/CAZy.clean.dmnd \\
            -q {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.threshold} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/blastp.sh'
        utils.write_cmd(cmd, shellname)
    

    def anno(self):
        # kegg数据库注释 
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8 \\
            --cazy_relation ~/pipeline/metagenomics/database/CAZy/cazy_relation.json \\
            --result {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8.anno 
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/anno.sh'
        utils.write_cmd(cmd, shellname)

    
    def abundance_stat(self):
        # 给cazy注释结果增加物种和丰度信息 
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_stat_abundance.py \\
            --cazy_anno {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8.anno  \\
            --abundance_table {self.projdir}/4.TaxAnnotation/NonRundant.tax_abundance.absolute.total.xls\\
            --sample_file  ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt\\
            --result_suffix {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8.anno.cazy
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/CAZy/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)


    def gene_stat(self):
        # 统计不同cazy水平的基因计数
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/cazy_stat_gene.py \\
            --cazy_anno {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum \\
            --sample_file ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt \\
            --result_suffix {self.projdir}/5.FunctionAnnotation/CAZy/NonRundant.protein.m8.anno.cazy.gene_stat
        
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
        self.projdir = args['projdir']
        self.threshold = self.args['threshold']
          

    def blastp(self):
        utils.mkdirs(f'{self.projdir}/5.FunctionAnnotation/eggNOG')
        cmd = textwrap.dedent(f'''
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d ~/pipeline/metagenomics/database/eggNOG/eggnog-mapper/data/e5.proteomes.dmnd  \\
            -q {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.support_reads.fa \\
            -o {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8 \\
            --max-target-seqs 1 \\
            --evalue {self.threshold} \\
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/blastp.sh'
        utils.write_cmd(cmd, shellname)


    def anno(self):
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_anno.py \\
            --m8 {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8 \\
            --eggNOG_relation ~/pipeline/metagenomics/database/eggNOG/eggnog-mapper/data/eggNOG_config/eggNOG.json \\
            --result {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8.anno 
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/anno.sh'
        utils.write_cmd(cmd, shellname)


    def abundance_stat(self):
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_stat_abundance.py \\
            --eggNOG_anno {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8.anno \\
            --abundance_table {self.projdir}/4.TaxAnnotation/NonRundant.tax_abundance.absolute.total.xls \\
            --sample_file ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt \\
            --result_suffix {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8.anno.eggNOG
        ''')
        shellname = f'{self.projdir}/5.FunctionAnnotation/eggNOG/abundance_stat.sh'
        utils.write_cmd(cmd, shellname)
         

    def gene_stat(self):
        cmd = textwrap.dedent(f'''
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/FunctionAnnotation/bin/eggNOG_stat_gene.py  \\
            --eggNOG_anno {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8.anno \\
            --gene_table {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum \\
            --sample_file ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/sample_info.txt \\
            --result_suffix {self.projdir}/5.FunctionAnnotation/eggNOG/NonRundant.protein.m8.anno.eggNOG.gene_stat\\
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
    

