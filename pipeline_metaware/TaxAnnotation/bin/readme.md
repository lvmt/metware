# 构造指定微生物NR子库

**包含物种信息**
 - 细菌：2
 - 真菌：4751
 - 古菌：2157
 - 病毒：10239

https://bioinf.shenwei.me/taxonkit/tutorial/


**运行方式**
   ```
   # 1. 列出所有id下面的taxid
   taxonkit list --ids 2,4751,2157,10239 --indent "" >  micro.taxid.txt
   
   # 2. Retrieving target accessions
   # time: 4min
    pigz -dc prot.accession2taxid.gz \
        | csvtk grep -t -f taxid -P $id.taxid.txt \
        | csvtk cut -t -f accession.version,taxid \
        | sed 1d \
        > $id.acc2taxid.txt

    cut -f 1 $id.acc2taxid.txt > $id.acc.txt

    # 3. 提取对应accession序列
    seqtk subseq nr.gz $id.acc.txt > $id.na.fa


   ```


