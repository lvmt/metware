# 宏基因组流程从头开发

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
    - 单样本预测(metagenemark)
    - unmap reads预测
    - 基因预测结果处理
        - 重命名基因名称(cd-hit不接受重复名)
        - 剔除核酸长度小于100nt的预测结果
    - CD-HIT预测(基于protein水平,主要是快些)
        - 使用大内存, 多线程;不然慢到你怀疑人生
    - 基于protein的CD-HIT预测结果,获取核酸的CD-HIT结果
    - bowtie2比对(clean fq与聚类后的nucleotide.fa)
    - 统计丰度信息(得到的为相对丰度)
    - 相对丰度转换为绝对丰度（理论上没有绝对丰度这一说法，其实就是在相对丰度矩阵上，乘以放大系数）

## 4、TaxAnnotation 物种注释目录
### 4.1 构建NR子库文件
**构造方式参考<a href="https://bioinf.shenwei.me/taxonkit/tutorial/#making-nr-blastdb-for-specific-taxids">txonkit</a>软件描述即可**

本流程主要构造了细菌，真菌，古菌，病毒四种微生物

1、了解上述四种物种的taxid
```
细菌：2
古菌：2157
真菌：4751
病毒：10239
```
2、获取每个物种下属的全部子taxid
```
taxonkit list --ids 10239 --indent "" > virus.taxid.txt
```
3、基于taxid，获取蛋白的accession
```
# prot.accession2taxid.gz需要下载
zcat prot.accession2taxid.gz | csvtk -t grep -f taxid -P bacteria.taxid.txt | csvtk -t cut -f accession.version,taxid > bacteria.taxid.accession.txt
```
4、从nr.fasta文件里面，提取特定accession的序列
```
seqtk subseq nr.gz bacteria.taxid.accession.txt > bacteria.na.fa
```
5、diamond构建索引
```
diamond makedb --in nr.fa -d nr  # 输出结果：nr.dmnd
```
### 4.2 物种注释步骤
#### 注释软件
    - diamond
    - MEGAN
```
1、 diamond blastp 
2、MEGAN daa2rma rma2info 提取需要信息
3、通过库文件，给上述结果基于taxid增加物种注释信息
4、合并物种注释和丰度统计文件，并拆分为不同层级
```




## 5、FunctionAnnotation 功能注释目录

## 6、database
**存放常用的数据库文件**



