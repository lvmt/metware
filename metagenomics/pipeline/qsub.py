#!/usr/bin/env python3
# -*- coding:utf-8 -*- 
'''
@Author: lvmengting 
@Date: 2022-04-28 15:57:15
'''

# PBS投递系统-简版-只是投递-无冗余功能

dos = """ 
\033[1;32m                                                                           
         _______                   _____                    _____                    _____          
        /::\    \                 /\    \                  /\    \                  /\    \         
       /::::\    \               /::\    \                /::\____\                /::\    \        
      /::::::\    \             /::::\    \              /:::/    /               /::::\    \       
     /::::::::\    \           /::::::\    \            /:::/    /               /::::::\    \      
    /:::/~~\:::\    \         /:::/\:::\    \          /:::/    /               /:::/\:::\    \     
   /:::/    \:::\    \       /:::/__\:::\    \        /:::/    /               /:::/__\:::\    \    
  /:::/    / \:::\    \      \:::\   \:::\    \      /:::/    /               /::::\   \:::\    \   
 /:::/____/   \:::\____\   ___\:::\   \:::\    \    /:::/    /      _____    /::::::\   \:::\    \  
|:::|    |     |:::|    | /\   \:::\   \:::\    \  /:::/____/      /\    \  /:::/\:::\   \:::\ ___\ 
|:::|____|     |:::|____|/::\   \:::\   \:::\____\|:::|    /      /::\____\/:::/__\:::\   \:::|    |
 \:::\   _\___/:::/    / \:::\   \:::\   \::/    /|:::|____\     /:::/    /\:::\   \:::\  /:::|____|
  \:::\ |::| /:::/    /   \:::\   \:::\   \/____/  \:::\    \   /:::/    /  \:::\   \:::\/:::/    / 
   \:::\|::|/:::/    /     \:::\   \:::\    \       \:::\    \ /:::/    /    \:::\   \::::::/    /  
    \::::::::::/    /       \:::\   \:::\____\       \:::\    /:::/    /      \:::\   \::::/    /   
     \::::::::/    /         \:::\  /:::/    /        \:::\__/:::/    /        \:::\  /:::/    /    
      \::::::/    /           \:::\/:::/    /          \::::::::/    /          \:::\/:::/    /     
       \::::/____/             \::::::/    /            \::::::/    /            \::::::/    /      
        |::|    |               \::::/    /              \::::/    /              \::::/    /       
        |::|____|                \::/    /                \::/____/                \::/____/        
         ~~                       \/____/                  ~~                       ~~              
                                                                                                    
                  -- tool for simplest qsub                                  
\033[0m
"""


from argparse import ArgumentDefaultsHelpFormatter
import os  
import subprocess
import textwrap 



class QSUB:

    def __init__(self, args):
        self.args = args
        self.sh = os.path.abspath(self.args['sh'])
        self.shell_abshpath = os.path.abspath(self.sh) 
        self.shell_absdir = os.path.dirname(self.sh)
        self.mem = args['mem']
        self.threads = args['threads']


    def get_stdout(self):
        return self.args['o'] or f'{self.sh}.stdout'

    def get_stderr(self):
        return self.args['o'] or f'{self.sh}.stderr'

    def get_task_name(self):
        return os.path.basename(self.sh)

    def get_mailer(self):
        mailer = subprocess.getoutput("who i am| awk '{print $1}'")
        return f'{mailer}.metware.cn'

    def get_priority(self):
        return float(self.args['p']) or 10

    def get_queue(self):
        # 给予队列的闲置情况,进行任务投递
        if self.args['q']:
            return self.args['q']
        node, queue = subprocess.getoutput('pestat  | grep free  | awk \'{if($2~"middle" || $2 ~ "high") print $1,$2}\'|head -1').split(' ')
        queue = queue.split(',')[-1]
        return queue


    def start(self):
        stdout = self.get_stdout()
        stderr = self.get_stderr()
        task_name = self.get_task_name()
        mailer = self.get_mailer()
        priority = self.get_priority()
        queue = self.get_queue()

        cmd = textwrap.dedent(f""" 
        #!/bin/bash
        #PBS -N {task_name} 
        #PBS -l select=1:ncpus={self.threads}:mem={self.mem}
        #PBS -q {queue}
        #PBS -p {priority}
        #PBS -M {mailer}
        #PBS -o {stdout}
        #PBS -e {stderr} 

        cd {self.shell_absdir}

        """)
        #PBS -e {stderr}  
        pbs_shellname = f'{self.sh}.pbs'
        with open(self.sh, 'r') as fr, open(pbs_shellname, 'w') as fw:
            fw.write(cmd)
            sh_text = fr.readlines()
            fw.write(''.join(sh_text))
         
        os.system(f'nohup ssh -YT mgr qsub {pbs_shellname} &')




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=dos, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('sh', help='\033[1;32m待运行程序, one shell one job\033[0m')
    parser.add_argument('-a', help='程序等待时间')
    parser.add_argument('-C', help='以该字符开头的行qsub的命令选项')
    parser.add_argument('-e', help='标准错误信息定向路径')
    parser.add_argument('-o', help='标准输出重定向路径')
    parser.add_argument('-j', help='将标准输出与标准错误合并输出到一个文件join中')
    parser.add_argument('-N', help='任务名称')
    parser.add_argument('-M', help='邮件发送给谁')
    parser.add_argument('-p', help='定义任务优先级[-1024, 1023]', type=int, default=0)
    parser.add_argument('-q', help='队列信息')
    parser.add_argument('--mem', help='内存大小2G,default=2gb', default='2gb')
    parser.add_argument('--threads', help='线程数目,default=1', default=1)

    args = vars(parser.parse_args())

    QSUB(args).start()





