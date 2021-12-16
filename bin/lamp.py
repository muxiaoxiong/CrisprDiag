#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-09 19:41:42
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com
#03版需增加自由能的检验!!!!符合条件的seq只有163种!!!!不太可行，故先加入回文序列检测，引物末端重复碱基检测

"""
LAMP引物设计
形成3个引物对
F3 F2 F1 B1 B2 B3
GC含量与Tm值的关系

		   >=60    <=45
F1c/B1c   64-66   59-61
F2/B2     59-61   54-56
F3/B3     59-61   54-56
LF/LB     59-61   54-56

引物最末端6个碱基 自由能<=-4K
F3/B3 F2/B2 LF/LB 3‘末端G<=-4k
F1c/B1c 5’末端G<=-4k

GC含量在40%-60%

F2 5‘--B2 5’   距离在120-180bp  (120~160bp)
F2 5'--F1c 3'  距离在40-60bp
F2 5'--F3  3'  距离在0-20bp    (0~60bp)

F3-->B3 长度为300bp

内引物 FIP BIP
FIP--> F1c+F2
BIP--> B1c+B2
外引物 F3 B3

总引物长度300bp
F3--B3 长度15-25bp

不能全部枚举

A complementary sequence is defined as symmetric
sequences (for example CCCGGG and GAATTC) and
special sequences (for example, sequences containing the 
same nucleotide at the end such as CCGGGG and
AATTTT)

The 5’ end of F1c after amplification corresponds to the 3’ end of F1, so that stability is 
important.
"""
import primer3
# from tools import CountGC,Fasta_reverse,get_dict
import random
from bin import tools

can_len_GC=[18,19,17,20,16,21,15,22]
can_len_AT=[21,22,20,23,19,24,18,25]
loop_len=[50,49,51,48,52,47,53,46,54,45,55,44,56,43,57,42,58,51,59,40,60]

def work(sgRNA_dict,gene_file,outfile,m,tm1,tm2,tm3,tm4,num=10):
	print('正在设计LAMP引物')
	out=open(outfile,'a+')
	out.write('sgRNAID\tseq\tstrand\tstart\tend\tGC\t')
	for i in range(m+1):
		out.write('m%s\t'%i)
	out.write('F3\tB3\tFIP\tBIP\tLF\tLB\n')
	gene_dict=tools.get_dict(gene_file)
	for key,value in sgRNA_dict.items():
		geneseq=''.join(gene_dict[key])
		lenseq=len(geneseq)
		sgRNAlist=value
		if lenseq<300:
			print('%s序列长度小于300bp,无法设计lamp引物...'%key)
		else:
			for sgRNAinfo in sgRNAlist:
				start=sgRNAinfo[3]
				end=sgRNAinfo[4]
				# print(sgRNAinfo[2])
				if start<150 or lenseq-end<150:
					print('%s无法设计lamp引物——该片段位置太前或太后'%sgRNAinfo[0])
				else:
					tar_gc=tools.CountGC(geneseq)
					if tar_gc<0.45:
						tar=1
						can_len=can_len_AT
					else:
						tar=0
						can_len=can_len_GC
					F1=''
					B1=''
					F2=''
					B2=''
					F3=''
					B3=''
					candidateF1=[]
					candidateB1C=[]
					candidate=[]
					if sgRNAinfo[2]=="-":
						for i1 in can_len:#F1长度
							for j1 in range(20):  # ？？FIP是否可以与sgRNA区域重叠
								j1=random.randint(0,19)
								F1start=end-j1-i1
								F1end=F1start+i1
								F1seq=geneseq[F1start:F1end]
								F1_primer=tools.Fasta_reverse(F1seq)#F1引物需反向互补，则B1引物不需反向互补
								j1=random.randint(0,19)
								B1Cstart=end+j1
								B1Cend=B1Cstart+i1
								B1Cseq=geneseq[B1Cstart:B1Cend]
								B1_primer=B1Cseq
								if select_prime1(tar,F1_primer,tm1,tm2) == 1 and JudgeSelfComplementary(F1_primer,maxnumber=4) and Judge_repeat(F1_primer) and dg(F1_primer[:7])<-4:
										                 
									candidateF1.append([F1_primer,F1start,F1end])
								if select_prime1(tar,B1_primer,tm1,tm2) == 1 and JudgeSelfComplementary(B1_primer,maxnumber=4) and Judge_repeat(B1_primer) and dg(B1_primer[:7])<-4:
									                     
									candidateB1C.append([B1_primer,B1Cstart,B1Cend])
					if sgRNAinfo[2]=="+":
						for i1 in can_len:#F1长度
							for j1 in range(20):  # ？？FIP是否可以与sgRNA区域重叠
								j1=random.randint(0,19)
								F1start=start-j1-i1
								F1end=F1start+i1
								F1seq=geneseq[F1start:F1end]
								F1_primer=tools.Fasta_reverse(F1seq)#F1引物需反向互补，则B1引物不需反向互补
								j1=random.randint(0,19)
								B1Cstart=start+j1
								B1Cend=B1Cstart+i1
								B1Cseq=geneseq[B1Cstart:B1Cend]
								B1_primer=B1Cseq
								if select_prime1(tar,F1_primer,tm1,tm2) == 1 and JudgeSelfComplementary(F1_primer,maxnumber=4) and Judge_repeat(F1_primer) and dg(F1_primer[:7])<-4:
										                 
									candidateF1.append([F1_primer,F1start,F1end])
								if select_prime1(tar,B1_primer,tm1,tm2) == 1 and JudgeSelfComplementary(B1_primer,maxnumber=4) and Judge_repeat(B1_primer) and dg(B1_primer[:7])<-4:
									                     
									candidateB1C.append([B1_primer,B1Cstart,B1Cend])
						#确定F2 和B2
					if len(candidateF1)==0 or len(candidateB1C)==0:
						print('%s无法设计lamp引物——无符合条件的F1、B1'%sgRNAinfo[0])
					else:
							#确定F2
						count=0
						count_F1=0
						for i in range(len(candidateF1)):
							n=random.randint(0,len(candidateF1)-1)
							i=candidateF1[n]
							for j in range(len(candidateB1C)):
								n=random.randint(0,len(candidateB1C)-1)
								j=candidateB1C[n]
								F1seq=i[0]
								F1start=i[1]
								F1end=i[2]
								B1Cseq=j[0]
								B1Cstart=j[1]
								B1Cend=j[2]
								for i1 in can_len:
									for j1 in loop_len:
										F2start=F1start-j1
										F2end=F2start+i1
										F2seq=geneseq[F2start:F2end]
										F2_primer=F2seq
										if select_prime2(tar,F2_primer,tm3,tm4) == 1 and JudgeSelfComplementary(F2_primer,maxnumber=4) and Judge_repeat(F2_primer) and dg(F2_primer[-6:])<-4:
											for i2 in can_len:#F2长
												for j2 in loop_len:#指针向后移
													B2Cend=B1Cend+j2
													B2Cstart=B2Cend-i2
													B2Cseq=geneseq[B2Cstart:B2Cend]
													B2_primer=tools.Fasta_reverse(B2Cseq)
													if select_prime2(tar,B2_primer,tm3,tm4) == 1 and F2start-B2Cend<180 and JudgeSelfComplementary(B2_primer,maxnumber=4) and Judge_repeat(B2_primer) and dg(B2_primer[-6:])<-4:
														result=[F1seq,F1start,F1end,B1Cseq,B1Cstart,B1Cend,F2_primer,F2start,F2end,B2_primer,B2Cstart,B2Cend]
														candidate.append(result)
														count+=1
														break
													if count>=num:
														break
												if count>=num:
													break
										if count>=num:
											break
									if count>=num:
										break
								if count>=num:
									break
							if count>=num:
								break
						if len(candidate)==0:
							print('%s无法设计lamp引物——无符合条件的F2、B2'%sgRNAinfo[0])
						else:              
							for i in candidate:
								F1_primer=i[0]
								F2_primer=i[6]
								B1_primer=i[3]
								B2_primer=i[9]
								F1start=i[1]
								F2start=i[7]
								B1Cend=i[5]
								F2end=i[8]
								B2Cstart=i[10]
								B2Cend=i[11]
								flag=0
								for i3 in can_len:
									for j3 in range(20):
										j3=random.randint(0,19)
										F3_primer=""
										F3start=F2start-j3-i3
										F3end=F3start+i3
										F3seq=geneseq[F3start:F3end]									
										if select_prime2(tar,F3seq,tm3,tm4)==1 and JudgeSelfComplementary(F3seq,maxnumber=4) and Judge_repeat(F3seq) and dg(F3seq[-6:])<-4:
											F3_primer=F3seq
											flag=1
											break
									if flag !=0:
										break
								flag=0
								for i3 in can_len:
									for j3 in range(20):
										j3=random.randint(0,19)
										B3_primer=""
										B3Cstart=B2Cend+j3
										B3Cend=B3Cstart+i3
										B3Cseq=geneseq[B3Cstart:B3Cend]
										B3seq=tools.Fasta_reverse(B3Cseq)
										if select_prime2(tar,B3seq,tm3,tm4)==1 and JudgeSelfComplementary(B3seq,maxnumber=4) and Judge_repeat(B3seq) and dg(B3seq[-6:])<-4:
											B3_primer=B3seq
											flag=1
											break
									if flag!=0:
										break
								if B3_primer !='' and F3_primer !='':
									flag=0
									LB=""
									for i4 in can_len:
										LB=""
										if i4<B2Cstart-B1Cend:
											for j4 in range(B2Cstart-B1Cend-i4):
												j4=random.randint(0,B2Cstart-B1Cend-i4-1)
												LB=""
												LBstar=B1Cend+1+j4
												LBend=LBstar+i4
												LBseq=geneseq[LBstar:LBend]
												if select_prime2(tar,LBseq,tm3,tm4)==1 and JudgeSelfComplementary(LBseq,maxnumber=4) and Judge_repeat(LBseq):
													LB=LBseq
													flag=1
													break
											if flag==1:
												break
									flag=0
									LF=""
									for i4 in can_len:
										LF=""
										if i4<F1start-F2end:
											for j4 in range(F1start-F2end-i4):												
												j4=random.randint(0,F1start-F2end-i4-1)
												LF=""
												LFstar=F2end+1+j4
												LFend=LFstar+i4
												# print(LFstar-F2end)
												LFseq=tools.Fasta_reverse(geneseq[LFstar:LFend])
												if select_prime2(tar,LFseq,tm3,tm4)==1 and JudgeSelfComplementary(LFseq,maxnumber=4) and Judge_repeat(LFseq):

													LF=LFseq
													flag=1
													break
											if flag==1:
												break
									if LB !='' and LF !='':
										FIP=F1_primer+F2_primer
										BIP=B1_primer+B2_primer
										sgRNAinfo=[str(i) for i in sgRNAinfo]
										out.write('\t'.join(sgRNAinfo)+'\t')
										out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(F3_primer,B3_primer,FIP,BIP,LF,LB))
									else:
										FIP=F1_primer+F2_primer
										BIP=B1_primer+B2_primer
										sgRNAinfo=[str(i) for i in sgRNAinfo]
										out.write('\t'.join(sgRNAinfo)+'\t')
										out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(F3_primer,B3_primer,FIP,BIP,"*","*"))
								else:
									print('%s无法设计lamp引物——无符合条件的F3、B3'%sgRNAinfo[0])
	out.close()	                    
						
def select_prime1(tar,seq,t1,t2):
	gc=tools.CountGC(seq)
	if  tar==1:
		if 0.4<=gc<=0.65:
			tm=primer3.calcTm(seq)
			if t1<=tm<=t2:
				return 1
			else:
				return 0
		else:
			return 0
	else:
		if 0.4<=gc<=0.65:
			tm=primer3.calcTm(seq)
			if t1+5<=tm<=t2+5:
				return 1
			else:
				return 0
		else:
			return 0

def select_prime2(tar,seq,t3,t4):
	gc=tools.CountGC(seq)
	if  tar==1:
		if 0.4<=gc<=0.65:
			tm=primer3.calcTm(seq)
			if t3<=tm<=t4:
				return 1
			else:
				return 0
		else:
			return 0
	else:
		if 0.4<=gc<=0.65:
			tm=primer3.calcTm(seq)
			if t3+5<=tm<=t4+5:
				return 1
			else:
				return 0
		else:
			return 0

def dg(seq):
    dg=primer3.calcHomodimer(seq).dg
    return dg

def JudgeSelfComplementary(seq,maxnumber=4):#？？？是否需要连续，不连续是否也会形成二级结构
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

def Judge_repeat(seq):
	flag=True
	for i in ["AAAA","TTTT","CCCC","GGGG"]:
		if seq[-4:]== i:
			flag=False
	return flag