#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-07 10:44:08
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-07 10:44:08
'''

'''爬取CAZy需要的信息
'''


import requests
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
import re 
from collections import defaultdict
import time 


url = 'http://www.cazy.org/GT1_characterized.html'


class Spider:

    def __init__(self, args):
        self.urls = args['urls']
        self.result = args['result']


    def get_response(self, url):
        ua = UserAgent()
        header = ua.random 
        headers = {
            'User-Agent': header
        }
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
            return response.text
        else:
            exit('返回码错误')


    def get_soup(self, response):
        try:
            return BeautifulSoup(response, 'html.parser')
        except Exception as e:
            return BeautifulSoup(response, 'xml.parser')


    def get_pages(self, soup):
        return ['http://www.cazy.org/' + page.get('href') for page in soup.select('table span.pages a')]

    
    def parser_soup(self, soup):
        '''
        acc : [ec, name, decs]
        '''
        access_info = {}
        description = soup.select('table >tr:nth-child(1)>td.tdsum ')[0].text
        if description.startswith(('Deleted family!', 'No known')):  # 有些家族已经被delete了
            return 'Delete'
        contents = soup.find_all('tr', {'valign':{"top"}})
        for content in contents:
            tds = content.select('td')
            name = tds[0].text.strip()
            ec = '|'.join([ec.text for ec in tds[1].select('a')])
            accs = re.findall(r'([A-Z0-9_]+.\d*)<', str(tds[4]).split('<b>')[-1])
            for acc in accs:
                access_info[acc] = [ec, name, description]
        return access_info


    def start(self):
        next_urls = []
        # all_info = {}
        with open(self.urls, 'r') as fr, open(self.result, 'w', encoding='utf-8') as fw:
            for line in fr:
                url = line.strip('\n')
                response = self.get_response(url)
                soup = self.get_soup(response)
                info = self.parser_soup(soup)
                if info == 'Delete':
                    continue
                pages = self.get_pages(soup)
                next_urls.extend(pages)

                for _key in info:
                    ec, name, desc = info[_key]
                    fw.write('{url}\t{_key}\t{ec}\t{name}\t{desc}\n'.format(**locals()))
                    time.sleep(1)

            for url in next_urls:
                response = self.get_response(url)
                soup = self.get_soup(response)
                info = self.parser_soup(soup)
                for _key in info:
                    ec, name, desc = info[_key]
                    fw.write('{url}\t{_key}\t{ec}\t{name}\t{desc}\n'.format(**locals()))
                time.sleep(1)






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='hahhaha')
    parser.add_argument('--urls', help='存放url的文件')
    parser.add_argument('--result', help='输出文件目录')

    args = vars(parser.parse_args())
    Spider(args).start()


