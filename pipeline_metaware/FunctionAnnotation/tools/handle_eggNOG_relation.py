#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-15 14:09:18
'''

#制作eggNOG 注释库文件 


from collections import defaultdict
import json 


def get_og_anno(annofile):
    # 后续分析过程中，发现部分COG本人人工注释与网页注释结果无法衔接
    '''
    {
        'OG': [
                [taxid, _class, descript],
                [taxid, _class, descript]
                
                ]
    }
    '''
    og_anno_info = defaultdict(list)
    with open(annofile, 'r') as fr:
        for line in fr:
            linelist = line.strip('\n').split('\t')
            taxid, OG, _class, descr, *_ = linelist + ['NA']  # 有的OG没有desc
            if OG.startswith('KOG'):
                continue
            if not OG.startswith(('COG', 'arCOG')):
                OG = f'ENOG50{OG}'

            og_anno_info[OG].append([taxid, _class, descr])

    return og_anno_info


def main(args):
    '''
    {
        'access': [OG1, OG2]
    }
    OG1 = 'OG1|_class|descr'
    '''
    json_info = defaultdict(list)
    og_anno_info = get_og_anno(args['annofile'])

    for memfile in args['member_files']:
        with open(memfile, 'r') as fr:
            for line in fr:
                linelist = line.strip('\n').split('\t')
                OG, access = linelist
                access_list = access.split(',')
                if OG.startswith('KOG'):
                    continue
                if not OG.startswith(('COG', 'arCOG')):
                    OG = f'ENOG50{OG}'

                for acc in access_list:
                    og_anno_list = og_anno_info[OG]
                    OG_taxid, OG_class, OG_desc = og_anno_list[0] if og_anno_list[0][0] == '1' else og_anno_list[-1]
                    json_info[acc].append('|'.join((OG, OG_taxid,OG_class, OG_desc)))

    with open(args['out_json'], 'w') as fw:
        json.dump(json_info, fw, indent=4)




    


if __name__ == '__main__':
    # 我们只构造细菌，古菌，真菌，病毒，真核生物的config文件
    import argparse
    parser = argparse.ArgumentParser(description='制作eggNOG库文件')
    parser.add_argument('--annofile', help='OG 描述信息')
    parser.add_argument('--member_files', nargs='+', help='各个物种的access与OG关系')
    parser.add_argument('--out_json', help='json配置文件')

    args = vars(parser.parse_args())
    main(args)




    