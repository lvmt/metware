trimmomatic \
	PE  \
	-threads 4 \
	HC1_L1_1.fq.gz \
	HC1_L1_2.fq.gz \
	HC1_L1_qc_1.fq.gz \
	HC1_L1_1_se \
	HC1_L1_qc_2.fq.gz \
	HC1_L1_2_se \
	ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
	LEADING:2 \
	TRAILING:2 SLIDINGWINDOW:4:2 \
	MINLEN:25



# 去除低质量序列