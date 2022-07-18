#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-23 10:36:35
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-23 10:36:35
'''

# 物种注释过程中的画图函数
from ctypes import util
import textwrap


class TaxonPlot:

    def __init__(self, args):
        self.args = args 
        self.projdir = self.args['projdir']
        self.metastat_compare = self.args['metastat_compare']
        self.random_forest_compare = self.args['random_forest_compare']
        
    
    @staticmethod
    def write_cmd(cmd, shellname):
        with open(shellname, 'w') as fw:
            fw.write(cmd)


    def plotPCA(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/plotPCA.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -s {self.projdir}/group.txt \\
                -o {self.projdir}/4.TaxAnnotation/Plot/PCA &&\\
            
            echo "end time: " `date`
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/pca_{_class}.sh'
            self.write_cmd(cmd, shellname)
            
    
    def plotHeatmap(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/plotHeatmap.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                --scale row \\
                -a both \\
                -n both \\
                --outdir {self.projdir}/4.TaxAnnotation/Plot/heatmap \\
                --prefix taxonomy.{_class}.heatmap  &&\\

            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/plotHeatmap.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -r TRUE \\
                -g {self.projdir}/group.txt \\
                --scale row \\
                -a both \\
                -n both \\
                --outdir {self.projdir}/4.TaxAnnotation/Plot/heatmap \\
                --prefix taxonomy.{_class}.heatmap  &&\\

            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/plotHeatmap.R \\
                -i {self.projdir}/4.TaxAnnotation/GeneStat/NonRundant.tax_gene.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                -a both \\
                -n both \\
                --cluster row \\
                --outdir {self.projdir}/4.TaxAnnotation/Plot/GeneNums.BetweenSamples.heatmap/{_class} \\
                --prefix {_class}.genenum.heatmap  &&\\

            echo "end time: "  `date`
            ''')
            shellname = shellname = f'{self.projdir}/4.TaxAnnotation/Plot/plotheatmap_{_class}.sh'
            self.write_cmd(cmd, shellname)
         

    def anosim(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/Anosim.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                -o {self.projdir}/4.TaxAnnotation/Plot/Anosim  &&\\
            
            echo "end time: " `date`
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/anosim_{_class}.sh'
            self.write_cmd(cmd, shellname)


    def krona(self):
        cmd = textwrap.dedent(f'''
        echo "start time: " `date`
        /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/krona.R \\
            -i {self.projdir}/4.TaxAnnotation/Absolute/NonRundant.tax_abundance.absolute.species.xls \\
            -g {self.projdir}/group.txt \\
            -k /share/software/apps/anaconda3/envs/microbiome/bin/ktImportText \\
            -o {self.projdir}/4.TaxAnnotation/Plot/Krona  &&\\
                
        echo "end time: " `date`
        ''')
        shellname = f'{self.projdir}/4.TaxAnnotation/Plot/krona.sh'
        self.write_cmd(cmd, shellname)
         

    def calc_nmds(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            export PATH=$PATH:/share/software/apps/anaconda3/envs/microbiome/bin &&\\
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/calc_NMDS.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                -o {self.projdir}/4.TaxAnnotation/Plot/NMDS  &&\\
            echo "end time: " `date`
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/nmds_{_class}.sh'
            self.write_cmd(cmd, shellname)
         

    def upgma_stacked_barplot(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/UPGMA_stacked_barplot.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                -o {self.projdir}/4.TaxAnnotation/Plot/UPGMA \\
                -d bray &&\\
            echo "end time: " `date`
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/upgma_{_class}.sh'
            self.write_cmd(cmd, shellname)
         

    def lefse(self):
        cmd = textwrap.dedent(f'''
        /share/software/apps/anaconda3/envs/microbiome/bin/python3 /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/LEfSe_input_prepare.py \\
            -i {self.projdir}/4.TaxAnnotation/Relative \\
            -g {self.projdir}/group.txt \\
            -o {self.projdir}/4.TaxAnnotation/Plot/LEfSe \\
            -c {self.projdir}/4.TaxAnnotation/04.TaxAnnotation.yaml &&\\

        /share/software/apps/anaconda3/envs/microbiome/bin/python3 /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/LEfSe_difference.py \\
            -i {self.projdir}/4.TaxAnnotation/Relative \\
            -g {self.projdir}/group.txt \\
            -o {self.projdir}/4.TaxAnnotation/Plot/LEfSe \\
            -c {self.projdir}/4.TaxAnnotation/04.TaxAnnotation.yaml &&\\

        ''')
        shellname = f'{self.projdir}/4.TaxAnnotation/Plot/lefse.sh'
        self.write_cmd(cmd, shellname)


    def metastat(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            for compare in self.metastat_compare:
                cmd = textwrap.dedent(f'''
                echo "start time: " `date`
                /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/species_metastat.R \\
                    -i {self.projdir}/4.TaxAnnotation/Absolute/NonRundant.tax_abundance.absolute.{_class}.xls \\
                    -g {self.projdir}/group.txt \\
                    -o {self.projdir}/4.TaxAnnotation/Plot/Metastats/{_class} \\
                    --diff {compare} \\
                    -m /share/software/apps/anaconda3/envs/microbiome/bin/mothur  &&\\
                echo "end time: " `date`
                ''')
                shellname = f'{self.projdir}/4.TaxAnnotation/Plot/metastat_{_class}_{compare}.sh'
                self.write_cmd(cmd, shellname)


    def random_forest(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            for compare in self.random_forest_compare:
                cmd = textwrap.dedent(f'''
                echo "start time: " `date`
                /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/random_forest.R  \\
                    --input {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                    --group {self.projdir}/group.txt \\
                    --outdir {self.projdir}/4.TaxAnnotation/Plot/random_forest/{_class} \\
                    --diff {compare}  &&\\
                echo "end time: " `date`
                ''')
                shellname = f'{self.projdir}/4.TaxAnnotation/Plot/random_forest_{_class}_{compare}.sh'
                self.write_cmd(cmd, shellname)
            

    def calc_pcoa(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            echo "start time: " `date`
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/calc_PCoA.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -g {self.projdir}/group.txt \\
                -d bray \\
                -o {self.projdir}/4.TaxAnnotation/Plot/PCoA &&\\
            echo "end time: " `date`
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/pcoa_{_class}.sh'
            self.write_cmd(cmd, shellname)
         

    def stack_barplot(self):
        for _class in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
            cmd = textwrap.dedent(f'''
            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/stackedBarplot.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -o {self.projdir}/4.TaxAnnotation/Plot/top \\
                --prefix {_class}.barplot \\
                -t 10  &&\\

            /share/software/apps/anaconda3/envs/R4.1.2/bin/Rscript /share/Bioinfo/tmp/tangliu/research/metagenome/bin/04.TaxAnnotation/bin/stackedBarplot.R \\
                -i {self.projdir}/4.TaxAnnotation/Relative/NonRundant.tax_abundance.relative.{_class}.xls \\
                -o {self.projdir}/4.TaxAnnotation/Plot/top_group \\
                --prefix {_class}.barplot \\
                -t 10  &&\\
            ''')
            shellname = f'{self.projdir}/4.TaxAnnotation/Plot/top_stack_barplot_{_class}.sh'
            self.write_cmd(cmd, shellname)
          

    def start(self):  # sourcery skip: use-dict-items, use-join
        plot_info = {
            'plotPCA': self.plotPCA,
            'plotHeatmap': self.plotHeatmap,
            'anosim': self.anosim,
            'krona': self.krona,
            'calc_nmds': self.calc_nmds,
            'upgma_stacked_barplot': self.upgma_stacked_barplot,
            'lefse': self.lefse,
            'metastat': self.metastat,
            'random_forest': self.random_forest,
            'calc_pcoa': self.calc_pcoa,
            'stack_barplot': self.stack_barplot
        }

        for _plot in plot_info:
            plot_info[_plot]()