#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-11-06 15:37:26
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

"""
PCR+sgRNA

sgRNA 设计
"""

import os
import sys
import gc
import platform
import pandas as pd
from itertools import islice

import regex as re

from bin import tools
from bin import pp
from bin import sgRNA
from bin import rpa
from bin import pcr
from bin import lamp

#判定键是否存在 用in dict速度更快
#循环字典 只循环一个速度更快

def main():
    gene_file,genome_file,outflag,pamseq,re_gene,re_genome,rv_gene,rv_genome,gcmax,gcmin,pamori,m,pcrb,rpab,lampb,tm1,tm2,tm3,tm4=pp.Input()
    #建立索引
    if(platform.system()=='Windows'):
        print('Windows系统')
    elif(platform.system()=='Linux'):
        print('Linux系统')
    print('--step1--sgRNA----')
    # #首先进行sgRNA 的设计 --->先设计
    sgRNAfile=outflag+'.sgRNA.txt'
    if os.path.exists(sgRNAfile):
        os.remove(sgRNAfile)
        print('已删除重复文件%s'%sgRNAfile)
    sgRNA_dict=sgRNA.work(re_gene,rv_gene,gene_file,gcmax,gcmin,pamori,pamseq,re_genome,rv_genome,genome_file,m)
    #sgRNA输出
    sgRNA.out(sgRNA_dict,sgRNAfile,m)

    #设计PCR引物
    if pcrb:
        pcrfile=outflag+'.pcr.txt'
        if os.path.exists(pcrfile):
            os.remove(pcrfile)
            print('已删除重复文件%s'%pcrfile)
        pcr.work(sgRNA_dict,gene_file,pcrfile,m)

    #设计RPA引物
    if rpab:
        rpafile=outflag+'.rpa.txt'
        if os.path.exists(rpafile):
            os.remove(rpafile)
            print('已删除重复文件%s'%rpafile)
        rpa.work(sgRNA_dict,gene_file,rpafile,m)

    #设计lamp
    if lampb:
        print('--step2--lamp----')
        lampfile=outflag+'.lamp.txt'
        if os.path.exists(lampfile):
            os.remove(lampfile)
            print('已删除重复文件%s'%lampfile)
        lamp.work(sgRNA_dict,gene_file,lampfile,m,tm1,tm2,tm3,tm4)



if __name__ == '__main__':
    main()
