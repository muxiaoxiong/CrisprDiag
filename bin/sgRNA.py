#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-09 15:33:35
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

'''
sgRNA部分的操作
保留活性判断的部分
'''
import regex as re
from bin import tools

def work(re_gene,rv_gene,gene_file,gcmax,gcmin,pamori,pamseq,re_genome,rv_genome,genome_file,m=3):
    global index_dict
    index_dict={}
    sgRNA_dict=get(re_gene,rv_gene,gene_file,gcmax,gcmin,pamori,pamseq)
    index_dict=index(re_genome,rv_genome,genome_file)
    print('正在计算脱靶...')
    if m==0:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m0(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    elif m==1:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m1(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    elif m==2:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m2(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    elif m==3:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m3(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    elif m==4:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m4(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    elif m==5:
        for gene,sgRNAlist in sgRNA_dict.items():
            for sgRNA in sgRNAlist:
                result=m5(sgRNA[1],pamseq,pamori)
                sgRNA.extend(result)
    return sgRNA_dict

def out(sgRNA_dict,sgRNAfile,m):
    out=open(sgRNAfile,'a+')
    out.write('gene\tsgRNAid\tseq\tstrand\tstart\tend\tgc\t')
    for i in range(m+1):
        out.write('m'+str(i)+'\t')
    out.write('\n')
    for key,value in sgRNA_dict.items():
        for i in value:
            out.write(key+'\t')
            for j in i:
                out.write(str(j)+'\t')
            out.write('\n')
    out.close()

def get(re_gene,rv_gene,gene_file,gcmax,gcmin,pamori,pamseq):
    #先设计sgRNA 只设计
    print('正在扫描基因文件序列...')
    gene_dict=tools.get_dict(gene_file)
    if pamori==3:
        print('正在进行sgRNA设计...')
        return pamori_3(re_gene,rv_gene,gene_dict,pamseq,gcmax,gcmin)
    elif pamori==5:
        print('正在进行sgRNA设计...')
        return pamori_5(re_gene,rv_gene,gene_dict,pamseq,gcmax,gcmin)

def pamori_3(re_gene,rv_gene,gene_dict,pam,gcmax,gcmin):
    re_obj=re.compile(re_gene)
    rv_obj=re.compile(rv_gene)
    out_dict={}
    for key,value in gene_dict.items():
        value=''.join(value)
        sgRNA_list_f=re_obj.finditer(value,overlapped=True)     #正义链
        sgRNA_list_s=rv_obj.finditer(value,overlapped=True)     #反义链
        count=0
        out_dict[key]=[]
        for sgRNA in sgRNA_list_f:
            seq=sgRNA.group()
            nseq=seq[:-len(pam)]
            gc=tools.CountGC(nseq)
            if gc<gcmin or gc>gcmax:
                pass
            else:
                count+=1
                name=key+'_s_'+str(count)
                strand='+'
                start=sgRNA.start()
                end=sgRNA.end()
                sgRNA_info=[name,seq,strand,start,end,gc]
                out_dict[key].append(sgRNA_info)
        count=0
        for sgRNA in sgRNA_list_s:
            seq=tools.Fasta_reverse(sgRNA.group())
            nseq=seq[:-len(pam)]
            gc=tools.CountGC(nseq)
            if gc<gcmin or gc>gcmax:
                pass
            else:
                count+=1
                name=key+'_a_'+str(count)
                strand='-'
                start=sgRNA.start()
                end=sgRNA.end()
                sgRNA_info=[name,seq,strand,start,end,gc]
                out_dict[key].append(sgRNA_info)
    return out_dict

def pamori_5(re_gene,rv_gene,gene_dict,pam,gcmax,gcmin):
    re_obj=re.compile(re_gene)
    rv_obj=re.compile(rv_gene)
    out_dict={}
    for key,value in gene_dict.items():
        out_dict[key]=[]
        value=''.join(value)
        sgRNA_list_f=re_obj.finditer(value,overlapped=True)     #正义链
        sgRNA_list_s=rv_obj.finditer(value,overlapped=True)     #反义链
        count=0
        for sgRNA in sgRNA_list_f:
            seq=sgRNA.group()
            nseq=seq[len(pam):]
            gc=tools.CountGC(nseq)
            if gc<gcmin or gc>gcmax:
                pass
            else:
                count+=1
                name=key+'_s_'+str(count)
                strand='+'
                start=sgRNA.start()
                end=sgRNA.end()
                sgRNA_info=[name,seq,strand,start,end,gc]
                out_dict[key].append(sgRNA_info)
        count=0
        for sgRNA in sgRNA_list_s:
            seq=tools.Fasta_reverse(sgRNA.group())
            nseq=seq[len(pam):]
            gc=tools.CountGC(nseq)
            if gc<gcmin or gc>gcmax:
                pass
            else:
                count+=1
                name=key+'_a_'+str(count)
                strand='-'
                start=sgRNA.start()
                end=sgRNA.end()
                sgRNA_info=[name,seq,strand,start,end,gc]
                out_dict[key].append(sgRNA_info)
    return out_dict

def index(re_genome,rv_genome,genome_file):
    print('正在对基因组建立索引...')
    genome_dict=tools.get_dict(genome_file)
    index_dict={}
    re_obj=re.compile(re_genome)
    rv_obj=re.compile(rv_genome)
    for value in genome_dict.values():
        value=''.join(value)
        sgRNA_list_f=re_obj.findall(value,overlapped=True)
        sgRNA_list_s=rv_obj.findall(value,overlapped=True)
        for sgRNA in sgRNA_list_f:
            if sgRNA not in index_dict:
                index_dict[sgRNA]=1
            else:
                index_dict[sgRNA]+=1
        for sgRNA in sgRNA_list_s:
            sgRNA=tools.Fasta_reverse(sgRNA)
            if sgRNA not in index_dict:
                index_dict[sgRNA]=1
            else:
                index_dict[sgRNA]+=1
    return index_dict

def m0(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0]
    if seq in index_dict:
        result[0]=index_dict[seq]
    return result

def m1(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0,0]
    seqlist=list(seq)
    maxlen=len(seqlist)
    if seq in index_dict:
        result[0]=index_dict[seq]
    for i1 in range(maxlen):
        temp1=['A','G','C','T']
        temp1.remove(seqlist[i1])
        for j1 in temp1:
            tempseq=seqlist.copy()
            tempseq[i1]=j1
            tempseq=''.join(tempseq)
            if tempseq in index_dict:
                result[1]+=index_dict[tempseq]
    return result

def m2(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0,0,0]
    seqlist=list(seq)
    maxlen=len(seqlist)
    if seq in index_dict:
        result[0]=index_dict[seq]
    for i1 in range(maxlen):
        temp1=['A','G','C','T']
        temp1.remove(seqlist[i1])
        for j1 in temp1:
            tempseq=seqlist.copy()
            tempseq[i1]=j1
            tempseq=''.join(tempseq)
            if tempseq in index_dict:
                result[1]+=index_dict[tempseq]
    for i1 in range(0,maxlen-1):
        for i2 in range(i1+1,maxlen):
            temp1=['A','G','C','T']
            temp1.remove(seqlist[i1])
            temp2=['A','G','C','T']
            temp2.remove(seqlist[i2])
            for j1 in temp1:
                for j2 in temp2:
                    tempseq=seqlist.copy()
                    tempseq[i1]=j1
                    tempseq[i2]=j2
                    tempseq=''.join(tempseq)
                    if tempseq in index_dict:
                        result[2]+=index_dict[tempseq]
    return result

def m3(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0,0,0,0]
    seqlist=list(seq)
    maxlen=len(seqlist)
    if seq in index_dict:
        result[0]=index_dict[seq]
    for i1 in range(maxlen):
        temp1=['A','G','C','T']
        temp1.remove(seqlist[i1])
        for j1 in temp1:
            tempseq=seqlist.copy()
            tempseq[i1]=j1
            tempseq=''.join(tempseq)
            if tempseq in index_dict:
                result[1]+=index_dict[tempseq]
    for i1 in range(0,maxlen-1):
        for i2 in range(i1+1,maxlen):
            temp1=['A','G','C','T']
            temp1.remove(seqlist[i1])
            temp2=['A','G','C','T']
            temp2.remove(seqlist[i2])
            for j1 in temp1:
                for j2 in temp2:
                    tempseq=seqlist.copy()
                    tempseq[i1]=j1
                    tempseq[i2]=j2
                    tempseq=''.join(tempseq)
                    if tempseq in index_dict:
                        result[2]+=index_dict[tempseq]
    for i1 in range(0,maxlen-3):
        for i2 in range(i1,maxlen-3):
            for i3 in range(i2,maxlen-1):
                for i4 in range(i2,maxlen):
                    temp1=['A','G','C','T']
                    temp1.remove(seqlist[i1])
                    temp2=['A','G','C','T']
                    temp2.remove(seqlist[i2])
                    temp3=['A','G','C','T']
                    temp3.remove(seqlist[i3])
                    temp4=['A','G','C','T']
                    temp4.remove(seqlist[i4])
                    for j1 in temp1:
                        for j2 in temp2:
                            for j3 in temp3:
                                for j4 in temp4:
                                    tempseq=seqlist.copy()
                                    tempseq[i1]=j1
                                    tempseq[i2]=j2
                                    tempseq[i3]=j3
                                    tempseq[i4]=j4
                                    tempseq=''.join(tempseq)
                                    if tempseq in index_dict:
                                        result[3]+=index_dict[tempseq]
    return result

def m4(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0,0,0,0,0]
    seqlist=list(seq)
    maxlen=len(seqlist)
    if seq in index_dict:
        result[0]=index_dict[seq]
    for i1 in range(maxlen):
        temp1=['A','G','C','T']
        temp1.remove(seqlist[i1])
        for j1 in temp1:
            tempseq=seqlist.copy()
            tempseq[i1]=j1
            tempseq=''.join(tempseq)
            if tempseq in index_dict:
                result[1]+=index_dict[tempseq]
    for i1 in range(0,maxlen-1):
        for i2 in range(i1+1,maxlen):
            temp1=['A','G','C','T']
            temp1.remove(seqlist[i1])
            temp2=['A','G','C','T']
            temp2.remove(seqlist[i2])
            for j1 in temp1:
                for j2 in temp2:
                    tempseq=seqlist.copy()
                    tempseq[i1]=j1
                    tempseq[i2]=j2
                    tempseq=''.join(tempseq)
                    if tempseq in index_dict:
                        result[2]+=index_dict[tempseq]
    for i1 in range(0,maxlen-3):
        for i2 in range(i1,maxlen-3):
            for i3 in range(i2,maxlen-1):
                for i4 in range(i2,maxlen):
                    temp1=['A','G','C','T']
                    temp1.remove(seqlist[i1])
                    temp2=['A','G','C','T']
                    temp2.remove(seqlist[i2])
                    temp3=['A','G','C','T']
                    temp3.remove(seqlist[i3])
                    temp4=['A','G','C','T']
                    temp4.remove(seqlist[i4])
                    for j1 in temp1:
                        for j2 in temp2:
                            for j3 in temp3:
                                for j4 in temp4:
                                    tempseq=seqlist.copy()
                                    tempseq[i1]=j1
                                    tempseq[i2]=j2
                                    tempseq[i3]=j3
                                    tempseq[i4]=j4
                                    tempseq=''.join(tempseq)
                                    if tempseq in index_dict:
                                        result[3]+=index_dict[tempseq]
    return result

def m5(seq,pam,pamori):
    if pamori==3:
        seq=seq[:-len(pam)]
    elif pamori==5:
        seq=seq[len(pam):]
    result=[0,0,0,0,0,0]
    seqlist=list(seq)
    maxlen=len(seqlist)
    if seq in index_dict:
        result[0]=index_dict[seq]
    for i1 in range(maxlen):
        temp1=['A','G','C','T']
        temp1.remove(seqlist[i1])
        for j1 in temp1:
            tempseq=seqlist.copy()
            tempseq[i1]=j1
            tempseq=''.join(tempseq)
            if tempseq in index_dict:
                result[1]+=index_dict[tempseq]
    for i1 in range(0,maxlen-1):
        for i2 in range(i1+1,maxlen):
            temp1=['A','G','C','T']
            temp1.remove(seqlist[i1])
            temp2=['A','G','C','T']
            temp2.remove(seqlist[i2])
            for j1 in temp1:
                for j2 in temp2:
                    tempseq=seqlist.copy()
                    tempseq[i1]=j1
                    tempseq[i2]=j2
                    tempseq=''.join(tempseq)
                    if tempseq in index_dict:
                        result[2]+=index_dict[tempseq]
    for i1 in range(0,maxlen-3):
        for i2 in range(i1,maxlen-3):
            for i3 in range(i2,maxlen-1):
                for i4 in range(i2,maxlen):
                    temp1=['A','G','C','T']
                    temp1.remove(seqlist[i1])
                    temp2=['A','G','C','T']
                    temp2.remove(seqlist[i2])
                    temp3=['A','G','C','T']
                    temp3.remove(seqlist[i3])
                    temp4=['A','G','C','T']
                    temp4.remove(seqlist[i4])
                    for j1 in temp1:
                        for j2 in temp2:
                            for j3 in temp3:
                                for j4 in temp4:
                                    tempseq=seqlist.copy()
                                    tempseq[i1]=j1
                                    tempseq[i2]=j2
                                    tempseq[i3]=j3
                                    tempseq[i4]=j4
                                    tempseq=''.join(tempseq)
                                    if tempseq in index_dict:
                                        result[3]+=index_dict[tempseq]
    for i1 in range(0,maxlen-4):
        for i2 in range(i1,maxlen-3):
            for i3 in range(i2,maxlen-2):
                for i4 in range(i2,maxlen-1):
                    for i5 in range(i2,maxlen):
                        temp1=['A','G','C','T']
                        temp1.remove(seqlist[i1])
                        temp2=['A','G','C','T']
                        temp2.remove(seqlist[i2])
                        temp3=['A','G','C','T']
                        temp3.remove(seqlist[i3])
                        temp4=['A','G','C','T']
                        temp4.remove(seqlist[i4])
                        temp5=['A','G','C','T']
                        temp5.remove(seqlist[i4])
                        for j1 in temp1:
                            for j2 in temp2:
                                for j3 in temp3:
                                    for j4 in temp4:
                                        for j5 in temp5:
                                            tempseq=seqlist.copy()
                                            tempseq[i1]=j1
                                            tempseq[i2]=j2
                                            tempseq[i3]=j3
                                            tempseq[i4]=j4
                                            tempseq[i5]=j5
                                            tempseq=''.join(tempseq)
                                            if tempseq in index_dict:
                                                result[4]+=index_dict[tempseq]
    return result


