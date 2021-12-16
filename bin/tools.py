#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-09 15:26:58
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

'''
序列方向互补
计算GC含量
'''
import primer3  #pip install primer3-py
import re
from bin import TmDeltaG as cal
# SantLucia 热力学参数表

def Fasta_reverse(sequence):
    #将序列进行方向互补
    sequence=sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]

def CountGC(seq):
    seq = seq.upper() #也可是使用.lower()把大写转换成小写计算
    count_c = seq.count('C')
    count_g = seq.count('G')
    gc_content = (count_g + count_c) / len(seq)
    return gc_content

def JudgeSelfComplementary(seq,maxnumber=4):
    #连续碱基互补数
    count=0
    seq1=Fasta_reverse(seq)
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
    if 45<=Tm1<=80 and 45<=Tm2<=80:
        return True
    else:
        return False

def judgeA(seq1,seq2):
    #判断碱基最后一位是否为A
    seq2=Fasta_reverse(seq2)
    if seq1[-1]=='A' or seq2[-1]=='A':
        return False
    else:
        return True

def judgeGGGCCC(seq1,seq2):
    #判断是否有3个以上G 或者 3个以上C
    seq2=fastaReverse(seq2)
    if seq1[-3:]=='GGG' or seq1[-3:]=='CCC' or seq2[-3:]=='GGG' or seq2[-3:]=='CCC':
        return False
    else:
        return True

def get_dict(file):
    a_dict={}
    with open(file) as ff:
        for line in ff:
            line=line.strip()
            if line.startswith('>'):
                name=line.replace('>','')
                a_dict[name]=[]
            else:
                line=line.upper()
                a_dict[name].append(line)
    return a_dict
