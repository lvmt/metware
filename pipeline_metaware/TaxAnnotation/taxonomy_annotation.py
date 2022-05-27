#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-20 09:51:38
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-20 09:51:38
'''


# 物种注释 taxonomy_annotation 
'''
1. diamond 上一步骤非荣誉蛋白序列与NR数据库比对(细菌,真菌,古菌,病毒)

'''

import re 
import os
import sys
import textwrap
import pandas as pd  

try:
    from Utils import utils
except Exception as e:
    sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
    from Utils import utils



class Taxonomy:

    def __init__(self, args):
        self.args = args 
        self.sample_file = self.args['sample_file']
        self.projdir = self.args['projdir']
        self.analysis_list = args['analysis_list']
        self.threshold = args['threshold']
        self.fq_info = utils.get_fq_info(self.sample_file)


    def diamond(self):
        # 利用非冗余的蛋白序列进行比对
        # 在非冗余的条件下,进行了reads支持数的筛选
        utils.mkdirs(f'{self.projdir}/4.TaxAnnotation')
        cmd = textwrap.dedent(f'''
        # diamond NR数据库比对; 输出格式为daa, 向后兼容megan
        # 后续结果处理,严重怀疑诺禾没有使用megan进行LCA分析,而是通过accesion-taxid-taxonomy之间的关系进行转换的
        ~/pipeline/metagenomics/software/diamond \\
            blastp \\
            -d ~/pipeline/metagenomics/database/microbe_nr/microbe_nr.dmnd \\
            -q {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.protein.support_reads.fa \\
            -o {self.projdir}/4.TaxAnnotation/NonRundant.protein.daa \\
            --outfmt 100 \\
            --evalue {self.threshold} 
        ''')
        shellname = f'{self.projdir}/4.TaxAnnotation/diamond.sh'
        utils.write_cmd(cmd, shellname)


    def megan(self):
        cmd = textwrap.dedent(f'''
        # MEGAN 物种注释;megan已经整合了lca的功能
        ~/.conda/envs/python2_lmt/bin/daa2rma  \\
            -i {self.projdir}/4.TaxAnnotation/NonRundant.protein.daa \\
            -ms 50 \\
            -me 0.01 \\
            -top 50 \\
            -mdb ~/pipeline/metagenomics/database/megan-map-Feb2022.db \\
            -o {self.projdir}/4.TaxAnnotation/NonRundant.protein.rma &&\\

        ~/.conda/envs/python2_lmt/bin/rma2info \\
            -i {self.projdir}/4.TaxAnnotation/NonRundant.protein.rma \\
            -r2c Taxonomy \\
            > {self.projdir}/4.TaxAnnotation/NonRundant.protein.info &&\\

        # 比对最终结果添加物种层级注释 
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/TaxAnnotation/bin/addtaxonomy.py \\
            --lca {self.projdir}/4.TaxAnnotation/NonRundant.protein.info \\
            --taxonomy_database ~/pipeline/metagenomics/database/microbe_nr/taxonomy/microbe.taxid_taxonomy \\
            --result {self.projdir}/4.TaxAnnotation/NonRundant.protein.taxonomy 

        ''')
        shellname = f'{self.projdir}/4.TaxAnnotation/megan.sh'
        utils.write_cmd(cmd, shellname) 


    def combine_abundance_taxonomy(self):
        # 整合物种丰度和物种注释结果  
        cmd = textwrap.dedent(f'''
        ## 整合绝对丰度表格, 并拆分成7个层级的结果
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/TaxAnnotation/bin/merge_abundance_taxonomy.py \\
            --abundance_table {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum.absolute.xls \\
            --taxonomy_table {self.projdir}/4.TaxAnnotation/NonRundant.protein.taxonomy \\
            --result_suffix  {self.projdir}/4.TaxAnnotation/NonRundant.tax_abundance.absolute &&\\


        ## 整合相对丰度结果, 并拆分成7个层级的结果
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/TaxAnnotation/bin/merge_abundance_taxonomy.py \\
            --abundance_table {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum.relative.xls \\
            --taxonomy_table {self.projdir}/4.TaxAnnotation/NonRundant.protein.taxonomy \\
            --result_suffix  {self.projdir}/4.TaxAnnotation/NonRundant.tax_abundance.relative 

        ## 整合物种注释和基因数目的结果, 在单样本和全部样本的水平上,拆分为7个层级结果
        python3 ~/gitlab/meta_genomics/metagenomics/pipeline_metaware/TaxAnnotation/bin/merge_gene_taxonomy.py \\
            --gene_table  {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.combine.readsNum\\
            --taxonomy_tabel {self.projdir}/4.TaxAnnotation/NonRundant.protein.taxonomy \\
            --result_suffix {self.projdir}/4.TaxAnnotation/NonRundant.tax_gene

        ''')
        shellname = f'{self.projdir}/4.TaxAnnotation/merge_abundance_taxonomy.sh'
        utils.write_cmd(cmd, shellname)


    def plot(self):
        # 绘图分析参数
        pass 
        


    def start(self):
        analysis_list = self.analysis_list.split(';')
        if '4' in analysis_list:
            # self.diamond()
            # self.megan()
            self.combine_abundance_taxonomy()





if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='物种注释步骤')
    parser.add_argument('--projdir', help='project analysis absolute dirname')
    parser.add_argument('--sample_file', help='样本信息配置文件')
    parser.add_argument('--analysis_list', help='具体分析步骤,后续待扩展', type=str, default='4')
    # diamond 参数
    parser.add_argument('--threshold', help='evalue阈值', type=float, default=0.00001)
    args = vars(parser.parse_args())
    print(args)

    Taxonomy(args).start()


    




