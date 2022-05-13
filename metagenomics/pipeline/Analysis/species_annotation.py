#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-09 15:58:23
'''

# 物种注释：基于diamond
'''
1 diamond比对
2 megan处理
'''


import textwrap
from Utils import utils


class LineageAnnotation:

    def __init__(self, args):
        self.args = args
        self.projdir = self.args['projdir']
        self.threshold = self.args['threshold']  # 默认值0.001


    def diamond(self):
        cmd = textwrap.dedent(f'''
        # diamond NR数据库比对; 输出格式为daa, 向后兼容megan
         ~/pipeline/metagenomics/software/diamond \\
            blastx \\
            -d bacteria.nr.dmnd \\
            -q {self.projdir}/3.GenePrediction/Cluster/NonRundant.total.nucleotide.fa \\
            -o {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.daa \\
            --outfmt 100 \\
            --evalue {self.threshold} &&\\
        ''')
        shellname = f'{self.projdir}/4.GeneAnnotation/anno/diamond.sh'
        utils.write_cmd(cmd, shellname)

    
    def megan(self):
        cmd = textwrap.dedent('''
        # MEGAN 物种注释
        daa2rma  \\
            -i {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.daa \\
            -ms 50 \\
            -me 0.01 \\
            -top 50 \\
            -mdb ~/pipeline/metagenomics/database/megan-map-Feb2022.db \\
            -o {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.rma &&\\

        rma2info \\
            -i {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.rma \\
            -r2c Taxonomy -n true -p true -v \\
            > {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy 

        # 物种注释结果整理
        python3 handle_taxonomy_file.py \\
            {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy \\
            {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy.simple  \\

        ''')
        shellname = f'{self.projdir}/4.GeneAnnotation/anno/megan.sh'
        utils.write_cmd(cmd, shellname)


    def get_taxonomy_abundance_table(self):
        # 获取物种-丰度表格
        # 获取物种-基因数目表格
        cmd = textwrap.dedent(f'''
        # 生成物种注释表格
        # 7种层级表格等
        python3 get_taxonomy_gene_abundantable.py \\
            --abundance_table {self.projdir}/3.GenePrediction/Cluster/total.abundance.tables \\
            --taxonomy_anno {self.projdir}/4.GeneAnnotation/anno/NonRundant.total.nucleotide.taxonomy.simple \\
            --out-dir {self.projdir}/4.GeneAnnotation/table/
        ''')
        shellname = f'{self.projdir}/4.GeneAnnotation/table/taxonomy_table.sh'
        utils.write_cmd(cmd, shellname)


    def plot(self):
        # 各种各样的绘图
        cmd = textwrap.dedent('''
        # krona 
        python3 krona_plot.py \\
            --infile {self.projdir}/4.GeneAnnotation/table/taxonomy.geneNums.table \\
            --outdir {self.projdir}/4.GeneAnnotation/krona\

        # 物种丰度前10绘图
        python barplot_top10.py \\
            --abundance_table taxonomy.abundance.table.s \\
            --abundance_stat total.taxonomy.abundance.table.s\\
            --out_dir \\
            --out-prefix \\


        
        ''')


    def start(self):
        self.diamond()
        
