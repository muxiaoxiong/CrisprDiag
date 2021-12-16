#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-10 10:18:34
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

"""
引物设计原则
引物与非特异扩增序列的同源性不要超过70%或有连续8个互补碱基同源。
引物长度为15~30bp，一般为20~27mer
G+C含量一般为40%~60%  Tm=58-60度  GC=30-80%。
3'端最后5个碱基内不能有多于2个的G或C.
引物中四种碱基的分布最好是随机的，不要有聚嘌呤或聚嘧啶的存在。尤其3′端不应超过3个连续的G或C
引物自身不应存在互补序列
引物3'端不能选择A，最好选择T。
"""

#输入需要的数目 进行筛选 不要全部计算

# from tools import get_dict
from bin import tools
import primer3
import re

len_can=[20,21,22,23,24,25,19,26,18,27,17,28,16,29,15,30]
def work(sgRNA_dict,gene_file,pcrfile,num=5):
    print('正在进行PCR引物设计...')
    gene_dict=tools.get_dict(gene_file)
    out=open(pcrfile,'a+')
    out.write('up\tstart\tend\tgc\tTm\tdown\tstart\tend\tgc\tTm\tcount\n')
    for key,value in gene_dict.items():
        if key in sgRNA_dict:
            geneseq=''.join(value)
            lenseq=len(geneseq)
            candidate=[]
            count=0
            for flag1 in range(lenseq):
                for flag2 in range(lenseq,1,-1):
                    for i1 in len_can:
                        for i2 in len_can:
                            up=geneseq[flag1:flag1+i1]
                            if JudgeGC(up) and Judge3(up) and Judge3GC(up) and JudgeSelfComplementary(up):
                                down=geneseq[flag2-i2:flag2]
                                if down=='':
                                    pass
                                else:
                                    if JudgeGC(down) and Judge3(down) and Judge3GC(down) and JudgeSelfComplementary(down):
                                        if judgeTm(up,down):
                                            upstart=flag1
                                            upend=flag1+i1
                                            downstart=flag2-i2
                                            downend=flag2
                                            candidate.append([up,upstart,upend,down,downstart,downend])
                                            count+=1
                                            if count>=num:
                                                break
                        if count>=num:
                            break
                    if count>=num:
                        break
                if count>=num:
                    break
            for i in candidate:
                up=i[0]
                down=i[3]
                upstart=i[1]
                upend=i[2]
                downstart=i[4]
                downend=i[5]
                out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t'%(up,upstart,upend,tools.CountGC(up),primer3.calcTm(up),down,downstart,downend,tools.CountGC(down),primer3.calcTm(down)))
                count=0
                for sgRNAinfo in sgRNA_dict[key]:
                    start=sgRNAinfo[3]
                    end=sgRNAinfo[4]
                    sgRNAid=sgRNAinfo[0]
                    if start>upend and end <downstart:
                        count+=1
                out.write('%s\n'%count)
    out.close()



def JudgeGC(seq):
    #GC含量
    gc=tools.CountGC(seq)
    if gc >=0.4 and gc <= 0.6:
        return True
    else:
        return False

def Judge3(seq):
    m3=seq[-1]
    if m3=='A':
        return False
    else:
        return True

def Judge3GC(seq):
    gc3=seq[-5:]
    result=re.findall(r'[GC]',gc3)
    if len(result) >2:
        return False
    else:
        return True

def JudgeSelfComplementary(seq,maxnumber=4):
    #连续碱基互补数
    count=0
    seq1=tools.Fasta_reverse(seq)
    seq1List=list(seq1)
    seq2List=list(seq)
    for i in range(len(seq)):
        if seq1List[i]==seq2List[i]:
            count+=1
        else:
            count=0
        if count>maxnumber:
            return False
    return True

def judgeTm(seq1,seq2):
    #Tm值判断
    Tm1=primer3.calcTm(seq1)
    Tm2=primer3.calcTm(seq2)
    if 55<=Tm1<=80 and 55<=Tm2<=80:
        return True
    else:
        return False

