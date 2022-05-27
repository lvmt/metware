#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-05-27 15:45:17
@Last Modified by:   lvmengting
@Last Modified time: 2022-05-27 15:45:17
'''


import requests
from bs4 import BeautifulSoup
from fake_useragent import UserAgent
import time 

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
        return self.get_request()
        


if __name__ == '__main__':
    import sys 
    infile = sys.argv[1]

    with open(infile, 'r') as fr: # T00001 文件
        for line in fr:
            line = line.strip('\n').split('\t')
            kid = line[0]
            url = f'http://rest.kegg.jp/list/{kid}'
            print(url)
            text = Spider(url).start()
            with open(f'{kid}.gene', 'w') as fw:
                fw.write(text)
                time.sleep(2)






