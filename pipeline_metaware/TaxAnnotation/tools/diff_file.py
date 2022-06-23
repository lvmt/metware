#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-23 08:36:53
@Last Modified by:   lvmengting
@Last Modified time: 2022-06-23 08:36:53
'''

# 由于分析方法的更新，对比更新前后，结果文件的差异  




from collections import defaultdict


class Diff:

    def __init__(self, args):
        self.args = args 
        self.file1 = self.args['file1']
        self.file2 = self.args['file2']
        self.result_suffix = self.args['result_suffix']


    @staticmethod
    def file2dict(file):
        dict_info = defaultdict(dict)
        with open(file, 'r') as fr:
            header = fr.readline().strip('\n').split('\t')
            
            for line in fr:
                linelist = line.strip('\n').split('\t')
                dict_info[linelist[0]] = dict(zip(header[1:], map(float, linelist[1:])))
        # print(dict_info)
        return dict_info

    @staticmethod
    def diff(d1, d2):
        return next(('diff' for key in d1 if d1[key] != d2.get(key)), 'same')


    def start(self):
        dict_info1 = self.file2dict(self.file1)
        dict_info2 = self.file2dict(self.file2) 
        diff_file1 = []
        diff_file2 = []

        all_keys = list(dict_info1.keys()) + list(dict_info2.keys())

        for key in all_keys:
            if self.diff(dict_info1.get(key, {}), dict_info2.get(key, {})) == 'diff':
                diff_file1.append([key] + list(dict_info1[key].values()) if dict_info1.get(key) else [key, 'none'])
                diff_file2.append([key] + list(dict_info2[key].values()) if dict_info2.get(key) else [key, 'none'])

        if diff_file1:
            with open(f'{self.result_suffix}.diff.{self.file1}', 'w') as fw:
                for row in diff_file1:
                    fw.write('{}\n'.format('\t'.join(map(str, row))))

        if diff_file2:
            with open(f'{self.result_suffix}.diff.{self.file2}', 'w') as fw:
                for row in diff_file2:
                    fw.write('{}\n'.format('\t'.join(map(str, row))))

        if not diff_file1 and not diff_file2:
            print('\033[1;32msame file\033[0m')


           



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='文件比对')
    parser.add_argument('--file1')
    parser.add_argument('--file2')
    parser.add_argument('--result_suffix')

    args = vars(parser.parse_args())
    Diff(args).start()
