# 迈维代谢宏基因组开发流程

## 1、QC 质控分析目录
### 肠道宏基因组
    - bowtie2去除宿主
    - fastp质控
    - clean reads信息统计
    - 全部样本生成统计结果

### 环境样本宏基因组
    - fastp质控
    - clean reads信息统计
    - 全部样本生成统计结果


## 2、Assembly 组装分析目录
    - MEGAHIT 单样本组装
    - bowtie2获取unmap reads
    - unmap reads混合组装
    - QUAST评估组装结果

## 3、GenePrediction 基因预测分析目录

## 4、TaxAnnotation 物种注释目录

## 5、FunctionAnnotation 功能注释目录

## 6、database
**存放常用的数据库文件**



