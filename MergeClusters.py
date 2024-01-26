#! usr/bin/python3
## Merge clusters post-GIANA
## To solve an explicit bug in GIANA: by only considering the nearest neighbor,
## GIANA may randomly divide a bigger, connected graph into several disconnected graphs
## Jan 25, 2024

import sys, os, re, resource
from os import path
import numpy as np
from Bio.Align import substitution_matrices
import time
from time import gmtime, strftime
from operator import itemgetter
from itertools import chain
from random import shuffle
from optparse import OptionParser
from collections import Counter
blosum62=substitution_matrices.load('BLOSUM62')

thr_s=3.7

inputFile = sys.argv[1]  ## input file is the output of GIANA
try:
    outputFile = sys.argv[2]  ## output file, optional
except IndexError:
    outputFile=inputFile+'_merged.txt'

blosum62n={}
for kk in blosum62.keys():
    a1=kk[0]
    a2=kk[1]
    vv=blosum62[kk]
    if vv>4:
        vv=4
    blosum62n[(a1,a2)]=vv
    if a1 != a2:
        blosum62n[(a2,a1)]=vv

def SeqComparison(s1,s2,gap=-6):
    n=len(s1)
    CorList=[]
    score=0
    for kk in range(0,n):
        aa=s1[kk]
        bb=s2[kk]
        if aa in ['.','-','*'] or bb in ['.','-','*']:
            if aa!=bb:
                score += gap
            continue
        if aa==bb:
#            score += min(4,blosum62[(aa,aa)])
            score += blosum62n[(aa,aa)]
            continue
        KEY=(aa,bb)
#        if KEY not in blosum62:
#            KEY=(bb,aa)
#        if KEY not in blosum62:
#            raise "Non-standard amino acid coding!"
        score+=blosum62n[KEY]
    return score

def NHLocalAlignment(Seq1,Seq2,gap_thr=1,gap=-6):
    n1=len(Seq1)
    n2=len(Seq2)
    if n1<n2:
        Seq=Seq1
        Seq1=Seq2
        Seq2=Seq
        nn=n2-n1
    else:
        nn=n1-n2
    if nn>gap_thr:
        return -1
    SeqList1=[Seq1]
    SeqList2=InsertGap(Seq2,nn)
    alns=[]
    SCOREList=[]
    for s1 in SeqList1:
        for s2 in SeqList2:
                SCOREList.append(SeqComparison(s1,s2,gap))
    maxS=max(SCOREList)
    return maxS

def InsertGap(Seq,n):
    ## Insert n gaps to Seq; n<=2
    if n==0:
        return [Seq]
    ns=len(Seq)
    SeqList=[]
    if(n==1):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            SeqList.append(SeqNew)
    if(n==2):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            for jj in range(0,ns+2):
                SeqNew0=SeqNew[0:jj]+'-'+SeqNew[jj:]
                SeqList.append(SeqNew0)
    return SeqList

def falign(s1, s2, V1, V2 ,st,VScore={}, UseV=True, gapn=1, gap=-6):
    mid1=s1[st:-2]
    mid2=s2[st:-2]
    if UseV:
        if V2==V1:
            V_score=4
        else:
            Vkey=(V1,V2)
            if Vkey not in VScore:
                Vkey=(V2,V1)
            if Vkey not in VScore:
                #print("V gene not found!")
                return 0
            else:
                V_score=VScore[Vkey]/20.0
    else:
        V_score=4.0
    aln=NHLocalAlignment(mid1,mid2,gapn,gap)
    score=aln/float(max(len(mid1),len(mid2)))+V_score
    return score



h = open(inputFile)
alines=h.readlines()
headlines=alines[0:2]
alines=alines[2:]  ## skip headerlines

clusterDict = {}
for ll in alines:
    ww = ll.strip().split('\t')
    cdr3 = ww[0]
    Vgene= ww[2]
    VgeneFamily = re.sub('-.+','',Vgene)
    info = '\t'.join(ww[2:])
    CL = ww[1]
    if CL not in clusterDict:
        clusterDict[CL]=[[cdr3, VgeneFamily, Vgene, info]]
    else:
        clusterDict[CL].append([cdr3, VgeneFamily, Vgene, info])

## Collapse cluster dictionary by length and V gene
vf=open('./VgeneScores.txt')  ## Use tcrDist's Vgene 80-score calculation
VScore={}
while 1:
    line=vf.readline()
    if len(line)==0:
        break
    ww=line.strip().split('\t')
    VScore[(ww[0],ww[1])]=int(float(ww[2]))/20
    VScore[(ww[1],ww[0])]=int(float(ww[2]))/20

SeqDict = {}
           
for kk in clusterDict:
    vv=clusterDict[kk]
    LL = len(vv[0][0])
    VF = vv[0][1]
    Key = str(LL)+'_'+VF
    if Key not in SeqDict:
        SeqDict[Key]=[vv]
    else:
        SeqDict[Key].append(vv)

## Compare each cluster within a given length and Vgene family
def IfMerge(CL1, CL2):
    nC1=len(CL1)
    nC2=len(CL2)
    merge_flag=0
    for ii in range(nC1):
        cl1=CL1[ii]
        ss1=cl1[0]
        vg1=cl1[2]
        for jj in range(nC2):
            cl2=CL2[jj]
            ss2=cl2[0]
            vg2=cl2[2]
            SCORE=falign(ss1, ss2, vg1, vg2, st=3, VScore=VScore)
            SCORE = SCORE/2
            if SCORE<=thr_s/2:
                #print(SCORE)
                return 0
            if SCORE>=thr_s:
                merge_flag=1
                break
        if merge_flag==1:
            break
    return merge_flag

def MergeAllClusters(idx, CLs):
    vvnew=[]
    for ii in idx:
        vv = CLs[ii]
        tmp = ['='.join(x) for x in vv]
        vvnew += tmp
    vvnew=list(set(vvnew))
    vvnew = [x.split('=') for x in vvnew]
    return vvnew

SeqDictNew = {}
t1=time.time()
for KK in SeqDict:
    print('-----Merging group %s -----' %KK)
    VV = SeqDict[KK]
    nV = len(VV)
    tagV = [-1]*nV
    count=1
    for ii in range(nV):
#        if ii % 1000 ==0:
#            print(" Processed %d sequences" %ii)
        if tagV[ii]>0:
            continue
        vvi = VV[ii]
        if len(vvi)>=100 or len(vvi)<=3:
            continue
        for jj in range(nV):
            if jj<=ii:
                continue
            vvj = VV[jj]
            if len(vvj)>=100 or len(vvj)<=3:
                continue
            mFlag = IfMerge(vvi, vvj)
            if mFlag==1:
                if tagV[ii]== -1:
                    tagV[ii]=count
                    tagV[jj]=count
                    count +=1
                else:
                    tagV[jj]=tagV[ii]
    tagV = np.array(tagV)
    nMerge = np.where((tagV== -1))[0]
    VVnew=[]
    for ss in nMerge:
        VVnew.append(VV[ss])
    nMax=np.max(tagV)
    for ss in range(nMax):
        vss = np.where((tagV==ss))[0]
        vvnew = MergeAllClusters(vss, VV)
        VVnew.append(vvnew)
    SeqDictNew[KK]=VVnew
    print(' Successfully merged %d clusters' %(len(VV)-len(VVnew)))
    print('     Elapsed %f seconds' %(time.time()-t1))

g = open(outputFile, 'w')
g.write(headlines[0])
g.write(headlines[1])
CL=0
for KK in SeqDictNew:
    VV=SeqDictNew[KK]
    for vv in VV:
        CL+=1
        for vs in vv:
            vsnew=[vs[0],str(CL), vs[2], vs[3]]
            ll = '\t'.join(vsnew)+'\n'
            g.write(ll)
g.close()



  























