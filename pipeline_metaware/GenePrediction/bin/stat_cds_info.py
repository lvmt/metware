#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-17 15:07:39
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-17 15:07:39
'''

# 对预测得到的核苷酸序列信息进行统计
'''
名称    长度    GC  起始密码子  终止密码子  前3个碱基   前3个碱基(反向互补链)   后3个碱基    后3个碱基(反向互补链)
'''
import Bio
from Bio import SeqIO
from Bio.SeqUtils import GC



def main(args):  # sourcery skip: remove-redundant-slice-index
    start_codon_pool = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}
    stop_codon_pool = {'TAA', 'TAG', 'TGA'}

    with open(args['result'], 'w') as fw:
        fw.write('Gene\tLength\tGC(%)\tInitiation_codon\tStop_codon\tFirst_3_letter\tFirst_3_letter(reverse complementary strand)\tLast_3_letter\tLast_3_letter(reverse complementary strand)\n')
        for seq_record in SeqIO.parse(args['fa'], 'fasta'):
            name = seq_record.id
            _seq = seq_record.seq 
            reverse_complement_seq = _seq.reverse_complement()
            length = len(seq_record)
            gc_percent = round(GC(_seq), 2)

            first_3_letter = _seq[:3]
            first_3_letter_rev_com = reverse_complement_seq[:3]
            last_3_letter = _seq[-3:]
            last_3_letter_rev_com = reverse_complement_seq[-3:]

            tmp_start = {first_3_letter, first_3_letter_rev_com}
            tmp_stop = {last_3_letter, last_3_letter_rev_com}

            if tmp_start & start_codon_pool:
                start_flag = 'yes'
            else:
                start_flag = 'no'

            if tmp_stop & stop_codon_pool:
                stop_flag = 'yes'
            else:
                stop_flag = 'no'

            fw.write(f'{name}\t{length}\t{gc_percent}\t{start_flag}\t{stop_flag}\t{first_3_letter}\t{first_3_letter_rev_com}\t{last_3_letter}\t{last_3_letter_rev_com}\n')
                


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='统计下基本信息关于预测得到的基因')
    parser.add_argument('--fa', help='输入文件: 待统计的fa')
    parser.add_argument('--result', help='输出文件: 统计结果')

    args = vars(parser.parse_args())
    main(args)
