#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-26 15:43:17
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-26 15:43:17
'''


# kegg api 工具接口使用



import requests
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
import textwrap



ua = UserAgent()
header = ua.random 
headers = {
    'User-Agent': header
}
print('ua', ua)



class Spider:

    def __init__(self, url):
        self.url = url
    
    def get_request(self):
        response = requests.get(self.url, headers=headers)
        if response.status_code == 200:
            return response.text
        print('wrong')
        exit(response.status_code)
        


    def get_soup(self, response):
        try:
            soup = BeautifulSoup(response, 'xml.parser')
        except Exception as e:
            soup = BeautifulSoup(response, 'html.parser')
        return soup 


    def start(self):
        response = self.get_request()
        return self.get_soup(response)




class INFO:
    
    BASE_URL  = 'http://rest.kegg.jp/info'

    def get_parser(self, parser):
        parser.add_argument('--db', 
                            help='展示哪个库',
                            choices=['pathway', 'ko', 'genome', 'enzyme', 'disease', 'kegg'])

        parser.set_defaults(func=self.start)

    def start(self, **args):
        db = args['db']
        print(f'query site: {self.BASE_URL}/{db}')
        print(Spider(f'{self.BASE_URL}/{db}').get_request())




class LIST:

    BASE_URL = 'http://rest.kegg.jp/list'
    # http://rest.kegg.jp/list/<database>[/<option>]

    def get_parser(self, parser):
        parser.add_argument('--db', 
                            choices=['pathway', 'ko', 'genome', 'organism', 'enzyme', 'disease', 'kegg'],
                            default='')
        parser.add_argument('--spices', help='某个具体的物种,如hsa',default='')
        parser.add_argument('--kid', help='具体的kid编号', default='')
        parser.add_argument('--gene', help='gene name', default='')
        parser.set_defaults(func=self.start)


    def start(self, **args):
        url = '{self.BASE_URL}/{db}/{kid}/{spices}'.format(**locals(), **args )
        if args['gene']:
            url = '{url}:{gene}'.format(**locals(), **args)
        print(f'\033[1;32mquery site: {url}\033[0m')
        print(Spider(f'{url}').start())



class LINK:

    BASE_URL = 'http://rest.kegg.jp/link/'

    def get_parser(self, parser):
        parser.add_argument('--target_db', 
                            choices=['pathway', 'ko', 'genome', 'enzyme', 'disease', 'kegg'],
                            default='')
        parser.add_argument('--source_db', 
                            choices=['pathway', 'ko', 'genome', 'enzyme', 'disease', 'kegg'],
                            default='')        
        parser.add_argument('--spices', help='某个具体的物种,如hsa',default='')
        parser.add_argument('--kid', help='具体的kid编号', default='')
        parser.add_argument('--gene', help='gene name', default='')
        parser.set_defaults(func=self.start)



    def start(self, **args):
        pass  



       
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='\033[1;32mKEGG API tools\033[0m'
    )

    # http://rest.kegg.jp/<operation>/<argument>[/<argument2[/<argument3> ...]]

    subparser = parser.add_subparsers(
        title = 'sub-command',
        dest='subparser_name',
        metavar=''
    )
    
    #info 子命令
    info = INFO()
    info_parser = subparser.add_parser('info',
                                        formatter_class=argparse.RawTextHelpFormatter,
                                        help='展示统计信息')
                                    
    info.get_parser(info_parser)


    # list 子命令
    _list = LIST()
    list_dos = textwrap.dedent('''\033[33;1m
    http://rest.kegg.jp/list/<database>[/<option>]
    \033[32;1mExamples
    /list/pathway 	   	    returns the list of reference pathways
    /list/pathway/hsa               returns the list of human pathways
    /list/organism 	   	    returns the list of KEGG organisms with taxonomic classification
    /list/hsa 	   	            returns the entire list of human genes
    /list/T01001 	   	    same as above
    /list/hsa:10458+ece:Z5100 	    returns the list of a human gene and an E.coli O157 gene
    /list/cpd:C01290+gl:G00092 	    returns the list of a compound entry and a glycan entry
    /list/C01290+G00092 	    same as above
    \033[0m''')
    list_parser = subparser.add_parser('list',
                                        formatter_class=argparse.RawTextHelpFormatter,
                                        description  = list_dos,
                                        help='列出查询条目包含内容')
                                    
    _list.get_parser(list_parser)


    # link 两个数据库的连接关系: 例如: 一个KO对应哪些基因
    link = LINK()
    link_dos = textwrap.dedent('''
    /link/pathway/hsa 	   	KEGG pathways linked from each of the human genes
    /link/hsa/pathway 	   	human genes linked from each of the KEGG pathways
    /link/pathway/hsa:10458+ece:Z5100 	   	KEGG pathways linked from a human gene and an E. coli O157 gene
    /link/genes/K00500 	   	List of genes with the KO assignment of K00500
    /link/hsa/hsa00010 	   	List of human genes in pathway hsa00010
    /link/ko/map00010 or /link/ko/ko00010 	   	List of KO entries in pathway map00010 or ko00010
    /link/rn/map00010 or /link/rn/rn00010 	   	List of reaction entries in pathway map00010 or rn00010
    /link/ec/map00010 or /link/ec/ec00010 	   	List of EC number entries in pathway map00010 or ec00010
    /link/cpd/map00010 	   	List of compound entries in pathway map00010
    /link/drug/ndc 	   	NDC codes linked from KEGG DRUG entries
    /link/ndc/D00564 	   	NDC codes that correspond to given a KEGG DRUG entry
    ''')
    link_parser = subparser.add_parser('link', 
                                        formatter_class=argparse.RawTextHelpFormatter,
                                        description=link_dos,
                                        help='展示查询条目在目标库中对应的结果')
    link.get_parser(link_parser)



    args = parser.parse_args()
    args.func(**vars(args))
