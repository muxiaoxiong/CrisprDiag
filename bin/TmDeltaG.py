#!/usr/bin/env python
# -*- coding: utf-8 -*-


import math
from bin import ThermodynamicsParameters as TP

def Fasta_reverse(sequence):
    #将序列进行方向互补
    sequence=sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]

def calDeltaHS(qseq, sseq):
    """ Calculate deltaH and deltaS """

    #qseq = re.sub('[atcgn]+', '', qseq)
    #sseq = re.sub('[atcgn]+', '', sseq)

    init_begin = 'init%s%s' % (qseq[0], sseq[0])
    init_end = 'init%s%s' % (qseq[-1], sseq[-1])

    if init_begin in TP.dH_full and init_end in TP.dH_full:
        deltaH = TP.dH_full[init_begin] + TP.dH_full[init_end]
        deltaS = TP.dS_full[init_begin] + TP.dS_full[init_end]
    else:
        deltaH = 0
        deltaS = 0

    for i in range(len(qseq) - 1):
        dinuc = '%s%s' % (qseq[i : (i+2)], sseq[i : (i+2)])
        if dinuc in TP.dH_full and dinuc in TP.dS_full:
            deltaH += TP.dH_full[dinuc]
            deltaS += TP.dS_full[dinuc]

    return deltaH, deltaS

def calDeltaG(qseq, sseq, mono_conc=50, diva_conc=1.5, dntp_conc=0.25, deltaH=None, deltaS=None):
    """ Calculate the free Gibbs energy """

    mono_conc = float(mono_conc)
    diva_conc = float(diva_conc)
    dntp_conc = float(dntp_conc)

    if not (deltaH and deltaS):
        deltaH, deltaS = calDeltaHS(qseq, sseq)

    # Calculate the free Gibbs energy
    tao = 273.15 + 37 # Constant temperature tao in Kelvin

    # Many thanks for the anonymous referee who help me fix the bug in last version.
    mono_conc = mono_conc + divalent2monovalent(diva_conc, dntp_conc)
    mono_conc = mono_conc / 1000

    deltaS_adjust = deltaS + 0.368 * (len(sseq) - 1) * math.log(mono_conc, math.e)

    deltaG = (deltaH * 1000 - tao * deltaS_adjust) / 1000
    return deltaG

def calTm(qseq, sseq, mono_conc=50, diva_conc=1.5, oligo_conc=50, dntp_conc=0.25, deltaH=None, deltaS=None):
    """ Calculate Tm value of amplicon"""

    mono_conc = float(mono_conc)
    diva_conc = float(diva_conc)
    oligo_conc = float(oligo_conc)
    dntp_conc = float(dntp_conc)

    if not (deltaH and deltaS):
        deltaH, deltaS = calDeltaHS(qseq, sseq)

    deltaH = deltaH * 1000

    oligo_conc = oligo_conc / 1000000000

    # Many thanks for the anonymous referee who help me fix the bug in last version.
    mono_conc = mono_conc + divalent2monovalent(diva_conc, dntp_conc)
    mono_conc = mono_conc / 1000

    deltaS = deltaS + 0.368 * (len(qseq) - 1) * math.log(mono_conc, math.e)

    Tm = deltaH / (deltaS + 1.987 * math.log(oligo_conc / 4, math.e)) - 273.15

    return Tm

def divalent2monovalent(divalent, dntp):
    '''Divalent to monovalent'''
    if divalent==0:
        dntp = 0

    if divalent < 0 or dntp <0:
        print >> sys.stderr, 'Error: divalent2monovalent, please contact Wubin Qu <quwubin@gmail.com>.'
        exit()

    if divalent < dntp:
        # According to the theory, melting temperature doesn't depend on divalent cations
        divalent = dntp

    return 120 * (math.sqrt(divalent - dntp))


class Cal():
    def __init__(self, qseq, sseq, mono_conc=50, diva_conc=1.5, oligo_conc=50, dntp_conc=0.25):
        deltaH, deltaS = calDeltaHS(qseq, sseq)
        self.Tm = calTm(qseq, sseq, mono_conc=mono_conc, diva_conc=diva_conc, oligo_conc=oligo_conc, dntp_conc=dntp_conc,deltaH=deltaH, deltaS=deltaS)
        self.DeltaG = calDeltaG(qseq, sseq, mono_conc=mono_conc, diva_conc=diva_conc, dntp_conc=dntp_conc,deltaH=deltaH, deltaS=deltaS)

def main():
    qseq = 'TTTGCAGGAAGTT'
    sseq = Fasta_reverse(qseq)
    # qseq = 'GGACACTCTATGGGAAAGAGTGTCC'
    # sseq = 'GGACACTCTATGGGAAAGAGTGTCC'

    mono = 50
    diva = 1.5
    oligo = 50
    dntp = 0.25
    seq = Cal(qseq, sseq, mono_conc=50, diva_conc=1.5, oligo_conc=50, dntp_conc=0.25)
    print('Tm: ', seq.Tm)
    print('DeltaG: ', seq.DeltaG)

if __name__ == '__main__':
    main()
