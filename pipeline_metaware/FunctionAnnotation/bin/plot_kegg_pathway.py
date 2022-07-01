#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-06-29 13:52:38
'''


# 给通路上色中指定的元素上色
# pillow  

'''
通过conf文件, 获取每个元素的位置信息 
借助pillow，对指定的元素进行颜色绘制 
'''

import pandas as pd 
from  PIL import Image, ImageDraw 
import matplotlib as mpl  
import numpy as np
import re  
import subprocess
import os  



def mkdir(dirname):
    if not os.path.exists(dirname):
        os.system(f'mkdir -p {dirname}')


def draw_line(im, box, color, width=2):
    draw = ImageDraw.Draw(im)
    if len(box) == 4:
        draw.line(box, fill=color, width=width)
    for i in range(0, len(box) - 3, 4):
        draw.line(box[i:i+4], fill=color, width=width)


def draw_rectangle_outline(im, box, color):
    draw = ImageDraw.Draw(im)
    draw.rectangle(box, outline = color)


def draw_rectangle_fill(im, box, color):
    region = im.crop(box)
    region_data = region.load()
    for x in range(region.size[0]):
        for y in range(region.size[1]):
            if region_data[x, y] == (255, 255, 255, 255):  # 背景色
                region_data[x, y] = color
    im.paste(region, box)


def draw_circle(im, box, color):
    draw = ImageDraw.Draw(im)
    draw.ellipse(box, fill=color)


def get_color_list(max_count):
    cmap = mpl.cm.get_cmap('Reds', max_count + 100) 
    color_buf = cmap(np.linspace(0, 1, max_count + 100))  # 颜色往后顺延一点, 因为起始颜色过于淡
    color_buf = color_buf[50:]
    return [list(map(int, x * 255)) for x in color_buf]  # 



def draw_single_map(map_conf, map_png, color_list, map_df, out_png):
    # 颜色列表, 根据基因数目的多少, 获取对应索引的颜色
    KO_list = map_df['KO'].tolist()
    conf_df = pd.read_csv(map_conf, sep='\t', header=None)
    conf_df[1] = conf_df[1].apply(lambda x: x.split('?')[1])
    conf_df = conf_df.drop(1, axis=1).join(conf_df[1].str.split('+', expand=True).stack().reset_index(level=1, drop=True).rename(1))
    conf_df.rename(columns={0: 'shape', 2: 'info', 1: 'KO'}, inplace=True)
    merge_df = pd.merge(map_df, conf_df, on='KO', how='inner')
    im = Image.open(map_png)
    for row in merge_df.values:
        ko, ko_count, ec_count, shape, info = row
        box = re.findall(r'\d+', shape)
        box = list(map(int, box))
        color = tuple(color_list[ko_count])

        if 'rect' in shape:
            # box = tuple(box)
            draw_rectangle_fill(im, box, color)
        elif 'circle' in shape:
            x, y, radius = box 
            box = [x - radius, y - radius, x + radius, y + radius]
            draw_circle(im, box, color)
        elif 'line' in shape:
            width = int(box[-1])
            box = tuple(box[:-1])
            color = (0, 0, 255, 255)
            draw_line(im, box, color, width)

    im.save(out_png)


def main(args):
    mkdir(args['result_dir'])
    df = pd.read_csv(args['kegg_tax'], sep='\t')
    description_df = pd.read_csv(args['pathway_description'], sep='\t', index_col=0)
    map_list = df['KO_Pathway'].unique().tolist()
    df = df[['KO_Pathway', 'KO', 'KO_EC']]
    
    with open('{result_dir}/stat.xls'.format(**args), 'w') as fw:
        fw.write('Pathway_ID\tPathway_Level1\tPathway_Level2\tPathway_Level3\tIdentified_ECs\n')
        for _map in map_list:
            # 绘制单张图片
            map_suffix = _map.replace('map', '')
            if not os.path.exists('{kegg_database}/ko{map_suffix}.png'.format(**args, **locals())):
                # print(_map)
                continue
            subprocess.getoutput('cp {kegg_database}/ko{map_suffix}.png  {result_dir}/'.format(**args, **locals()))
            subprocess.getoutput('cp {kegg_database}/ko{map_suffix}.conf {result_dir}/'.format(**args, **locals()))

            png = '{result_dir}/ko{map_suffix}.png'.format(**args, **locals())
            conf = '{result_dir}/ko{map_suffix}.conf'.format(**args, **locals())
            out_png = '{result_dir}/map{map_suffix}.png'.format(**args, **locals())

            map_df = df[df['KO_Pathway'] == _map]
            map_df = map_df.groupby('KO').count().reset_index()
            max_count = max(map_df['KO_Pathway'])
            color_list = get_color_list(max_count)

            draw_single_map(conf, png, color_list, map_df, out_png)

            ecs = '|'.join(set('|'.join(set(df[df['KO_Pathway'] == _map]['KO'].tolist())).split('|')))
            description = '\t'.join(description_df.loc[_map, :].tolist())
            fw.write(f'{_map}\t{description}\t{ecs}\n')



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='绘制通路注释图, 颜色深浅表示基因数的多少')
    parser.add_argument('--kegg_tax', help='pathway水平的基因统计文件')
    parser.add_argument('--pathway_description', help='每个pathway的描述信息pathway.levels.txt.name.modify', default='pathway.levels.txt.name.modify')
    parser.add_argument('--kegg_database', help='本地存储各种pathway图片的文件路径', default='/lustrefs/Database/KEGG/latest/maps')
    parser.add_argument('--result_dir', help='输出目录')


    args = vars(parser.parse_args())
    main(args)

   