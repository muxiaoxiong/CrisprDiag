#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-09 19:41:31
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

"""
RPA引物设计要求
长度：30-35 合适的可扩展到40
5’端的3-5个核苷酸应当避免G  C在这里是有益的
3’端 GC有利于结合
GC含量在 0.7-0.3之间
扩增区间 100-200bp 300bp以下
"""
# from tools import CountGC,get_dict,Fasta_reverse
from bin import tools

def work(sgRNA_dict,gene_file,outfile,m,rpanum=5):
    print('正在设计rpa引物...')
    out=open(outfile,'a+')
    out.write('geneid\tsgRNAID\tseq\tstrand\tstart\tend\tGC\t')
    for i in range(m+1):
        out.write('m%s\t'%i)
    for i in range(rpanum):
        j=i+1
        out.write('up-%s\tup-start-%s\tup-end-%s\tdown-%s\tdown-start-%s\tdown-end-%s\t'%(j,j,j,j,j,j))
    out.write('\n')
    gene_dict=tools.get_dict(gene_file)
    for key,value in sgRNA_dict.items():
        geneseq=''.join(gene_dict[key])
        lenseq=len(geneseq)
        sgRNAlist=value
        for sgRNAinfo in sgRNAlist:
            start=sgRNAinfo[3]
            end=sgRNAinfo[4]
            candidate=[]
            count=0
            if start<30 or lenseq-end<30 or lenseq<150:
                print('%s无法设计rpa引物'%sgRNAinfo[0])
            else:
                # 确定上游引物开始的位置
                if start<200:
                    pp=35
                else:
                    pp=start-200
                for flag1 in range(30,start):
                    for j1 in range(30,36):
                        if flag1-j1<0:
                            pass
                        else:
                            up=geneseq[flag1-j1:flag1]
                            if Judge_GC(up)>0 and Judge_repeat(up)<=5:
                                if Judge_5(up)>=0:
                                    for flag2 in range(end+1,300-flag1):
                                        for j2 in range(30,36):
                                            if flag2+j2>lenseq:
                                                pass
                                            else:
                                                down=geneseq[flag2:flag2+j2]
                                                if Judge_5(down)>=0 and Judge_3(down)>=0 and Judge_repeat(down)<=4:
                                                    result=[up,flag1-j1,flag1,down,flag2,flag2+j2]
                                                    candidate.append(result)
                                                    count+=1
                                                    if count>=rpanum:
                                                        break
                                        if count>=rpanum:
                                            break
                        if count>=rpanum:
                            break
                    if count>=rpanum:
                        break
                sgRNAinfo=[str(i) for i in sgRNAinfo]
                out.write('%s\t%s\t'%(key,'\t'.join(sgRNAinfo)))
                for rpa in candidate:
                    rpa=[str(i) for i in rpa]
                    out.write('\t'.join(rpa)+'\t')
                out.write('\n')
    out.close()


def Judge_5(seq):
    if seq[0]=='C':
        return 5
    elif seq[0]=='A':
        return 0
    elif seq[0]=='T':
        return 0
    elif seq[0]=='G':
        return -10

def Judge_3(seq):
    if seq[-1]=='G':
        return 5
    elif seq[-1]=='C':
        return 5
    elif seq[-1]=='T':
        return 0
    elif seq[-1]=='A':
        return 0

def Judge_repeat(seq):
    i=1
    templen=1
    maxlen=1
    while i<len(seq):
        if seq[i]==seq[i-1]:
            i+=1
            templen+=1
        else:
            i+=1
            templen=1
        if templen>maxlen:
            maxlen=templen
    return maxlen

def Judge_GC(seq):
    gc=tools.CountGC(seq)
    if gc<0.6 or gc>0.4:
        return 5
    elif 0.6<gc<0.7 or 0.3<gc<0.4:
        return 0
    else:
        return -5
