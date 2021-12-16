#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2020-12-09 15:37:00
# @Author  : Muxiaoxiong
# @email   : xiongweinie@foxmail.com

'''
获得用户输入
'''
import sys
import getopt


re_dict={'A':'A','G':'G','C':'C','T':'T','N':'[\S]',
         'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]',
         'B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]'
}
rv_dict={'A':'T','T':'A','G':'C','C':'G','N':'[\S]',
         'R':'[TC]','Y':'[GA]','S':'[CG]','W':'[TA]','K':'[CA]','M':'[TG]',
         'B':'[GCA]','D':'[TCA]','H':'[TGA]','V':'[TGC]'
}

def Input():
	#读取输入参数
    gene_file=''
    genome_file=''
    outflag=''
    pamlen=20
    pamseq='NGG'
    pamori=3
    gcmax=0.8
    gcmin=0.2
    m=3
    tm1=47
    tm2=61
    tm3=47
    tm4=56
    pcr=''
    rpa=''
    lamp=''
    if len(sys.argv) == 1:  # 无参数时 直接显示使用信息
        usage()
        sys.exit()
    else:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:g:o:s:r:l:m:a123',['help','gcmax=','gcmin=','tm1=','tm2=','tm3=','tm4='])
        for opt, value in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt=='-i':
                gene_file=value
            elif opt=='-g':
                genome_file=value
            elif opt=='-o':
                outflag=value
            elif opt=='-s':
                pamseq=value
            elif opt=='-r':
                pamori=int(value)
            elif opt=='-l':
                pamlen=int(value)
            elif opt=='-m':
                m=int(value)
            elif opt=='-a':
                print('请输入正确的参数\n --amax     The maximum value of GC content              <int>           [default:80]\n --amin     The minimum value of GC content              <int>           [default:20]')
                sys.exit()
            elif opt=='-1':
                pcr=True
            elif opt=='-2':
                rpa=True
            elif opt=='-3':
                lamp=True
            elif opt=='--gcmax':
                gcmax=int(value)
                if gcmax>1:
                    gcmax=gcmax/100
            elif opt=='--gcmin':
                gcmin=int(value)
                if gcmin>1:
                    gcmin=gcmin/100
            elif opt=='--tm1':
                tm1=int(value)
            elif opt=='--tm2':
                tm2=int(value)
            elif opt=='--tm3':
                tm3=int(value)
            elif opt=='--tm4':
                tm4=int(value)
    if gene_file=='':
        print('-i 没有输入参数')
        sys.exit()
    elif genome_file=='':
        print('-g 没有输入参数')
        sys.exit()
    elif outflag=='':
        print('-o 没有输入参数')
        sys.exit()
    else:
        pamseq_re=[re_dict[i] for i in list(pamseq)]
        pamseq_rv=[rv_dict[i] for i in list(pamseq)[::-1]]
        if pamori==3:
            re_gene='.{%s}%s'%(pamlen,''.join(pamseq_re))
            rv_gene='%s.{%s}'%(''.join(pamseq_rv),pamlen)
            re_genome='.{%s}(?=%s)'%(pamlen,''.join(pamseq_re))
            rv_genome='(?<=%s).{%s}'%(''.join(pamseq_rv),pamlen)
        else:
            re_gene='%s.{%s}'%(''.join(pamseq_re),pamlen)
            rv_gene='.{%s}%s'%(pamlen,''.join(pamseq_rv))
            re_genome='(?<=%s).{%s}'%(''.join(pamseq_re),pamlen)
            rv_genome='.{%s}(?=%s)'%(pamlen,''.join(pamseq_rv))
        print('----Information----')
        print('Gene File:%s'%gene_file)
        print('Genome File:%s'%genome_file)
        print('Out Flag:%s'%outflag)
        print('Pam seq:%s'%pamseq)
        print('Pamori:%s'%pamori)
        print('GC max:%s'%gcmax)
        print('GC min:%s'%gcmin)
        print('Mismatch:%s'%m)
        return gene_file,genome_file,outflag,pamseq,re_gene,re_genome,rv_gene,rv_genome,gcmax,gcmin,pamori,m,pcr,rpa,lamp,tm1,tm2,tm3,tm4

def usage():
    help_info="""
Usage:
    python3 RAVI.py  -i <gene.fa> -g <genome.fa> -o <outputflag> [options]
Options:
    -h         Display help information
    -s         PAM sequence                                 <string>         [default:NGG]
    -r         pamori.enter 5 or 3                           <int>           [default:3]
    -l         Length of protospacer                         <int>           [default:20]
    -m         The num of mismatch                           <int>           [default:3]
    --amax     The maximum value of GC content              <int>           [default:80]
    --amin     The minimum value of GC content              <int>           [default:20]
    -tm1       XXXX
    -tm2       XXXX
    -tm3       XXXX
    -tm4       XXXX
    -1         PCR
    -2         RPA
    -3         LAMP



PAMcode:
    Code        Base         Code        Base         Code        Base
     A         Adenine        R          A or G        B          C or G or T
     C         Cytosine       Y          C or T        D          A or G or T
     G         Guanine        S          G or C        H          A or C or T
     T         Thymine        W          A or T        V          A or C or G
     N         any base       K          G or T
                              M          A or C
    """
    print(help_info)

