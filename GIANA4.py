#! usr/bin/python3
## Isometric distance alignment algorithm to process ~100M CDR3 sequences
## Ultra-fast alignment with Blosum62 scoring
## Key algorithm: convert each CDR3 with length L into a L-dimensional space with PCA encoding
## Use approximate nearest neighbor search to remove sequences without any close neighbor
## Apply pairwise comparison to the remaining CDR3s to identify clusters
## updated on May 27, 2020. Added V gene and SW alignment
## Remaining issues:
##  1. Different lengths will group together with larger thr_iso values (>11)
##  2. V gene and SW score matrix may fall into smaller cliques that require dps.
## updated on May 29, 2020. Change SW alignment into k-mer guided to accelerate
## updated on Jun 2, 2020. Add variable gene clustering unit to solve the 2nd remaining issue
## Future improvement:
##  1. rewrite in c++
##  2. find a better isometric representation of blosum62
##  3. find a solution to process gap

## updated on Jun 29, 2020. Collapse identical TCRs.

## July 16, 2020: Change name into GIANA: Geometric Isometry based ANtigen-specific tcr Alignment
## Aug 24, 2020: Add GPU option
## Sep 26, 2020: Find a bug in identical CDR3 handling when V genes are different

import sys, os, re, resource
from os import path
import numpy as np
#from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align import substitution_matrices
import time
from time import gmtime, strftime
from operator import itemgetter
from itertools import chain
from random import shuffle
from optparse import OptionParser
from collections import Counter
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import faiss
from query import *

blosum62=substitution_matrices.load('BLOSUM62')
AAstring='ACDEFGHIKLMNPQRSTVWY'
AAstringList=list(AAstring)
cur_dir=os.path.dirname(os.path.realpath(__file__))+'/'

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

bl62={'A':[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
      'R':[-1,4,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
      'N':[-2,0,4,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
      'D':[-2,-2,1,4,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
      'C':[0,-3,-3,-3,4,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
      'Q':[-1,1,0,0,-3,4,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
      'E':[-1,0,0,2,-4,2,4,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
      'G':[0,-2,0,-1,-3,-2,-2,4,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
      'H':[-2,0,1,-1,-3,0,0,-2,4,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
      'I':[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
      'L':[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
      'K':[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,4,-1,-3,-1,0,-1,-3,-2,-2],
      'M':[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,4,0,-2,-1,-1,-1,-1,1],
      'F':[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,4,-4,-2,-2,1,3,-1],
      'P':[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,4,-1,-1,-4,-3,-2],
      'S':[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
      'T':[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,4,-2,-2,0],
      'W':[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,4,2,-3],
      'Y':[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,4,-1],
      'V':[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]}

bl62c=np.array([np.array(x) for x in list(bl62.values())])
bl62c=4-bl62c

embedding=MDS(n_components=13, n_init=100, max_iter=1000, eps=0.00001, dissimilarity='precomputed')
X=embedding.fit_transform(bl62c) 

bl62np={}
vkk=list(bl62.keys())
for ii in range(20):
    kk=vkk[ii]
    bl62np[kk]=np.array(list(X[ii,])+[0]*17)

   
AAencodingDict={}
for ii in range(len(AAstringList)):
    aa=AAstringList[ii]
    CODE=[0]*(ii)+[1]+[0]*(20-ii)
    AAencodingDict[aa]=np.array(CODE)

Ndim=16  ## optimized for isometric embedding
n0=Ndim*6
#M0=np.concatenate((np.concatenate((ZERO,M1),axis=1),np.concatenate((M1, ZERO),axis=1)))
ZERO=np.zeros((Ndim,Ndim))
II=np.eye(Ndim)
M0=np.concatenate((np.concatenate((ZERO,ZERO, II),axis=1),np.concatenate((II, ZERO, ZERO),axis=1),np.concatenate((ZERO,II, ZERO),axis=1)))
## Construct 6-th order cyclic group
ZERO45=np.zeros((Ndim*3,Ndim*3))
M6=np.concatenate((np.concatenate((ZERO45,M0),axis=1),np.concatenate((M0, ZERO45),axis=1)))

X=np.array([[-0.31230882, -0.53572156, -0.01949946, -0.12211268, -0.70947917,
        -0.42211092,  0.02783931,  0.02637933, -0.41760305,  0.21809875,
         0.53532768,  0.04833016,  0.07877711,  0.50464914, -0.26972087,
        -0.52416842],
       [ 0.29672002,  0.29005364,  0.18176298, -0.05103382, -0.34686519,
         0.58024228, -0.49282931,  0.62304281, -0.09575202,  0.30115555,
         0.09913529,  0.1577466 , -0.94391939, -0.10505925,  0.05482389,
         0.38409897],
       [-0.42212537,  0.12225749,  0.16279646,  0.60099009,  0.19734216,
         0.42819919, -0.33562418,  0.17036334,  0.4234109 ,  0.46681561,
        -0.50347222, -0.37936876,  0.1494825 ,  0.32176759,  0.28584684,
         0.68469861],
       [ 0.18599294, -0.44017825, -0.4476952 ,  0.34340976,  0.44603553,
         0.40974629, -0.60045935, -0.09056728,  0.22147919, -0.33029418,
         0.55635594, -0.54149972,  0.05459062,  0.57334159, -0.06227118,
         0.65299872],
       [-0.19010428,  0.64418792, -0.85286762,  0.21380295,  0.37639516,
        -0.67753593,  0.38751609,  0.55746524,  0.01443766,  0.1776535 ,
         0.62853954, -0.15048523,  0.55100206, -0.21426656,  0.3644061 ,
        -0.0018255 ],
       [ 0.7350723 ,  0.10111267,  0.55640019, -0.18226966,  0.51658102,
        -0.19321508, -0.46599027, -0.02989911,  0.4036196 , -0.11978213,
        -0.29837524, -0.30232765, -0.36738065, -0.1379793 ,  0.04362871,
         0.33553714],
       [ 0.41134047,  0.13512443,  0.62492322, -0.10120261, -0.03093491,
         0.23751917, -0.68338694,  0.05124762,  0.41533821,  0.46669353,
         0.31467277, -0.02427587,  0.15361135,  0.70595112, -0.27952632,
         0.32408931],
       [-0.33041265, -0.43860065, -0.5509376 , -0.04380843, -0.35160935,
         0.25134855,  0.53409314,  0.54850824,  0.59490287,  0.32669345,
        -0.45355268, -0.56317041, -0.55416297,  0.18117841, -0.71600849,
        -0.08989825],
       [-0.40366849,  0.10978974,  0.0280101 , -0.46667987, -0.45607028,
         0.54114052, -0.77552923, -0.10720425,  0.55252091, -0.34397153,
        -0.59813694,  0.15567728,  0.03071009, -0.02176143,  0.34442719,
         0.14681541],
       [ 0.19280422,  0.35777863,  0.06139255,  0.20081699, -0.30546596,
        -0.56901549, -0.15290953, -0.31181573, -0.74523217,  0.22296016,
        -0.39143832, -0.16474685,  0.58064427, -0.77386654,  0.19713107,
        -0.49477418],
       [-0.16133903,  0.22112761, -0.53162136,  0.34764073, -0.08522381,
        -0.2510216 ,  0.04699411, -0.25702389, -0.8739765 , -0.24171728,
        -0.24370533,  0.42193635,  0.41056913, -0.60378211, -0.65756832,
         0.0845203 ],
       [-0.34792144,  0.18450939,  0.77038332,  0.63868511, -0.06221681,
         0.11930421,  0.04895523, -0.22463059, -0.03268844, -0.58941354,
         0.11640045,  0.32384901, -0.42952779,  0.58119471,  0.07288662,
         0.26669673],
       [ 0.01834555, -0.16367754,  0.34900298,  0.45087949,  0.47073855,
        -0.37377404,  0.0606911 ,  0.2455703 , -0.55182937, -0.20261009,
         0.28325423, -0.04741146,  0.30565238, -0.62090653,  0.17528413,
        -0.60434975],
       [-0.55464981,  0.50918784, -0.21371646, -0.63996967, -0.37656862,
         0.27852662,  0.3287838 , -0.56800869,  0.23260763, -0.20653106,
         0.63261439, -0.22666691,  0.00726302, -0.60125196,  0.07139961,
        -0.35086639],
       [ 0.94039731, -0.25999326,  0.43922549, -0.485738  , -0.20492235,
        -0.26005626,  0.68776626,  0.57826888, -0.05973995, -0.1193658 ,
        -0.12102433, -0.22091354,  0.43427913,  0.71447886,  0.32745991,
         0.03466398],
       [-0.13194625, -0.12262688,  0.18029209,  0.16555524,  0.39594125,
        -0.58110665,  0.16161717,  0.0839783 ,  0.0911945 ,  0.34546976,
        -0.29415349,  0.29891936, -0.60834721,  0.5943593 , -0.29473819,
         0.4864154 ],
       [ 0.40850093, -0.4638894 , -0.39732987, -0.01972861,  0.51189582,
         0.10176704,  0.37528519, -0.41479418, -0.1932531 ,  0.54732221,
        -0.11876511,  0.32843973, -0.259283  ,  0.59500132,  0.35168375,
        -0.21733727],
       [-0.50627723, -0.1973602 , -0.02339884, -0.66846048,  0.62696606,
         0.60049717,  0.69143364, -0.48053591,  0.17812208, -0.58481821,
        -0.23551415, -0.06229112,  0.20993116, -0.72485884,  0.34375662,
        -0.23539168],
       [-0.51388312, -0.2788953 ,  0.00859533, -0.5247195 , -0.18021544,
         0.28372911,  0.10791359,  0.13033494,  0.34294013, -0.70310089,
        -0.13245433,  0.48661081,  0.08451644, -0.69990992,  0.0408274 ,
        -0.47204888],
       [ 0.68546275,  0.22581365, -0.32571833,  0.34394298, -0.43232367,
        -0.5041842 ,  0.04784017, -0.53067936, -0.50049908,  0.36874221,
         0.22429186,  0.4616482 ,  0.11159174, -0.26827959, -0.39372848,
        -0.40987423]])

bl62np={}
vkk=list(bl62.keys())
for ii in range(20):
    kk=vkk[ii]
    bl62np[kk]=np.array(list(X[ii,])+[0]*Ndim*5)

def EncodingCDR3(s, M, n0):
    sL=list(s)
    x=np.array([0]*n0)
    for ii in range(len(sL)):
        x = np.dot(M, (x+bl62np[sL[ii]]))
    return x

def BuildLengthDict(seqs, sIDs, vGene=[], INFO=[]):
    LLs=[10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
    LengthD={}
    SeqD={}
    VgeneD={}
    InfoD={}
    AAs=set(list(AAencodingDict.keys()))
    NAs=len(AAencodingDict)
    cNAs=0
    for ii in range(len(seqs)):
        ID=sIDs[ii]
        ss=seqs[ii]
        ssAA=set(list(ss))
        TMP=list(ssAA | AAs)
        if len(TMP) > NAs:
            ## CDR3 containing non amino acid letter
            #print('Warning: CDR3: '+ss + ' contains non amino acid letter!')
            cNAs+=1
            continue
        if len(vGene)>0:
            vv=vGene[ii]
        if len(INFO)>0:
            info=INFO[ii]
        L=len(ss)
        if L not in LLs:
            continue
        if L not in LengthD:
            LengthD[L]=[ID]
            SeqD[L]=[ss]
            if len(vGene)>0:
                VgeneD[L]=[vv]
            if len(INFO)>0:
                InfoD[L]=[info]
        else:
            LengthD[L].append(ID)
            SeqD[L].append(ss)
            if len(vGene)>0:
                VgeneD[L].append(vv)
            if len(INFO)>0:
                InfoD[L].append(info)
    if cNAs>0:
        print("Warning: Skipped %d sequences with non AA letter!" %(cNAs))
    return LengthD, VgeneD, InfoD, SeqD

def CollapseUnique(LD, VD, ID, SD):
    kks=LD.keys()
    LDu={}
    VDu={}
    IDu={}
    SDu={}
    for kk in kks:
        vvL=list(LD[kk])
        if len(VD)>0:
            vvV=list(VD[kk])
        else:
            vvV=['TRBV2-1*01']*len(vvL)
        vvI=list(ID[kk])
        vvS=list(SD[kk])
        zz=zip(vvL, vvS, vvV, vvI)
        zzs=sorted(zz, key = lambda x: (x[1], x[2]))
        nz=len(zzs)
        pointer_pre=0
        pointer_cur=1
        s_pre=zzs[pointer_pre][1]
        v_pre=zzs[pointer_pre][2]
        uS=[s_pre]
        uV=[v_pre]
        uI=[[zzs[pointer_pre][3]]]
        while pointer_cur < nz:
            s_cur=zzs[pointer_cur][1]
            v_cur=zzs[pointer_cur][2]
            if s_cur == s_pre and v_cur == v_pre:
                uI[len(uI)-1].append(zzs[pointer_cur][3])
                pointer_cur += 1
                continue
            else:
                uS.append(s_cur)
                uV.append(v_cur)
                uI.append([zzs[pointer_cur][3]])
                s_pre=s_cur
                v_pre=v_cur
                pointer_pre=pointer_cur
                pointer_cur += 1
        uL=[x for x in range(len(uS))]
        LDu[kk]=uL
        SDu[kk]=uS
        if len(VD)>0:
            VDu[kk]=uV
        IDu[kk]=uI
    return LDu, VDu, IDu, SDu


class CDR3:
    def __init__(self, s, sID, KS, st, ed):
        ## initialize with an input sequence
        ## s: input CDR3 sequence starting from C and ending with the first F in FGXG
        ## sID: unique identifier (increasing integers) given to each CDR3 sequence. Even identical CDR3s should have distinct sIDs
        ## KS: Kmer size
        ## st: the first 0:(st-1) amino acids will not be included in K-merization
        ## ed: the last L-ed amino acids will be skipped
        self.s=s
        self.ID=sID
        L=len(s)
        self.L=L
        sub_s=s[st: (L-ed)]
        Ls=len(sub_s)
        Kmer=[sub_s[x:(x+KS)] for x in range(0,Ls-KS+1)]
        self.Kmer=Kmer

class KmerSet:
    ## Kmer set for fast read searching based on mismatch-allowed Kmer index
    def __init__(self, Seqs, sIDs, KS, st, ed):
        ## initialize with a list of CDR3s, parse each CDR3 into Kmers, build Kmer-sID dictionary
        ## Seqs and sIDs must have the same length
        if len(Seqs) != len(sIDs):
            raise "Sequence and ID lists have different length. Please check input."
        KmerDict={}
        N=len(Seqs)
        self.N=N
        CDR3Dict={}
        LLs=[]
        for ii in range(0,N):
            s=Seqs[ii]
            sID=sIDs[ii]
            cc=CDR3(s,sID,KS,st,ed)
            CDR3Dict[cc.ID]=cc.Kmer
            KK=cc.Kmer
            LLs.append(cc.L)
            for kk in KK:
                if kk not in KmerDict:
                    KmerDict[kk]=[sID]
                else:
                    KmerDict[kk].append(sID)
        self.KD=KmerDict
        self.KS=KS
        self.CD=CDR3Dict
        self.LL=LLs
    def FindKmerNeighbor(self,kk):
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            kk_pre=[kk[0:jj]]*20
            kk_suf=[kk[(jj+1):KS]]*20
            kkn=list(zip(kk_pre,AAstringList,kk_suf))
            KS_n1+=[''.join(list(x)) for x in kkn]
        return KS_n1
    def FindKmerNeighbor2(self,kk):
        ## KS>=6, allowing 2 mismatches. CDR3 length must be >= 10
        KS=self.KS
        KS_n1=[]
        for jj in range(KS):
            for ii in range(KS):
                if ii<=jj:
                    continue
                kk_pre=[kk[0:jj]]*20
                kk_mid=[kk[(jj+1):ii]]*20
                kk_suf=[kk[(ii+1):KS]]*400
                kkn=list(zip(kk_pre,AAstringList,kk_mid))
                kkn=[''.join(list(x)) for x in kkn]
                kkn=[[x]*20 for x in kkn]
                kkn=list(chain(*kkn))
                kkn2=list(zip(kkn, AAstringList*20, kk_suf))
                kkn2=[''.join(list(x)) for x in kkn2]
                KS_n1+=kkn2
        return KS_n1
    def KmerIndex(self):
        ## For each K-mer, find its nearest neighbor with 1 character mismatch
        KKs=list(self.KD.keys())
        KS=self.KS
        KKs_set=set(KKs)
        Skk='_'.join(KKs)
        KI_Dict={}
        for kk in KKs:
##            kk_neighbor=[]
##            for jj in range(KS):
##                kk_pre=kk[0:jj]
##                kk_suf=kk[(jj+1):KS]
##                pat=kk_pre+'['+AAstring+']{1}'+kk_suf
##                p=re.compile(pat)
##                mm=[m.group() for m in p.finditer(Skk)]
##                kk_neighbor+=mm
            KS_n=set(self.FindKmerNeighbor(kk))
            kk_neighbor = KS_n & KKs_set
            KI_Dict[kk]=list(kk_neighbor)
        return KI_Dict
    def updateKD(self, KI):
        ## group sequences sharing motifs with 1-2 mismatches
        KD=self.KD
        KDnew={}
        for kk in KD:
            kkm=KI[kk]
            vvL=itemgetter(*kkm)(KD)
            if isinstance(vvL[0],list):
                vvL=list(chain(*vvL))
            KDnew[kk]=vvL
        return KDnew

def GenerateMotifGraph(mD,seqs,seqID):
    SeqShareGraph={}
    mDL={}
    for kk in mD:
        vv=mD[kk]
        LL=[]
        for v in vv:
            LL.append(len(seqs[v]))
        mDL[kk]=LL
    for kk in mD:
        vv=mD[kk]
        LL=mDL[kk]
        nv=len(vv)
        for ii in range(0,nv):
            id_1=vv[ii]
            L1=LL[ii]
            for jj in range(ii,nv):
                if jj==ii:
                    continue
                id_2=vv[jj]
                L2=LL[jj]
                if L2 != L1:
                    continue
                if id_1 not in SeqShareGraph:
                    SeqShareGraph[id_1]=[id_2]
                elif id_2 not in SeqShareGraph[id_1]:
                    SeqShareGraph[id_1].append(id_2)
                if id_2 not in SeqShareGraph:
                    SeqShareGraph[id_2]=[id_1]
                elif id_1 not in SeqShareGraph[id_2]:
                    SeqShareGraph[id_2].append(id_1)
    return SeqShareGraph

def generateSSG(Kset, CDR3s, k_thr=2):
    KD=Kset.KD
    KI=Kset.KmerIndex()
    KDnew=Kset.updateKD(KI)
    CD=Kset.CD
    LL=np.array(Kset.LL)
    SSG={}
    for kk in CD:
        vv=itemgetter(*CD[kk])(KDnew)
        if isinstance(vv[0],list):
            vv=list(chain(*vv))
        vv1=[]
        c=Counter(vv)
        for k in c:
            if c[k]>=k_thr:
                vv1.append(k)
        vv1=np.array(vv1)
        if len(vv1)==0:
            continue
        cdr3=CDR3s[kk]
        L0=len(cdr3)
        idx=np.where(LL[vv1]==L0)[0]
        if len(idx)==0:
            continue
        vvs=list(vv1[idx])
        vvs.remove(kk)
        if len(vvs)>0:
            SSG[kk]=vvs
    return SSG

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

def UpdateSSG(SSG, seqs, Vgenes, Vscore={}, UseV=True, gap=-6, gapn=1, cutoff=7.5):
    SSGnew={}
    count=0
    t1=time.time()
    N=len(list(chain(*list(SSG.values()))))
#    print("Number of pairs to be processed: %d" %N)
    for kk in SSG:
        s1=seqs[kk]
        V1=Vgenes[kk]
        VV=SSG[kk]
        for vv in VV:
            s2=seqs[vv]
            V2=Vgenes[vv]
            score=falign(s1, s2, V1, V2, st=3, VScore=Vscore, UseV=UseV, gap=-6, gapn=1)
            count+=1
            if count % 1000000 ==0:
                t2=time.time()
#                print("Processed %d pairs. Elapsed time %f" %(count, t2-t1))
            if score>=cutoff:
                if kk not in SSGnew:
                    SSGnew[kk]=[vv]
                else:
                    SSGnew[kk].append(vv)
    return SSGnew

def dfs(graph, start):
    '''
    Non-resursive depth first search
    '''
    visited = set()
    stack = [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(set(graph[vertex]) - visited)
    
    return visited

def IdentifyMotifCluster(SSG):
    ## Input SeqShareGraph dictionary representation of sparse matrix
    POS=set(SSG.keys())
    NP=len(POS)
    ClusterList=[]
    tmpL=set(chain(*ClusterList))
    count=0
    while 1:
            xx=POS ^ tmpL
            if len(xx)==0:
                break
            for ii in xx:
#            STACK=LoadComm([],ii)
                STACK=dfs(SSG,ii)
                tmpL = tmpL | STACK
                ClusterList.append(list(STACK))
#                tmpL=set(chain(*ClusterList))
                count+=1
                if count % 200 ==0:
                    print ("    Solved %d clusters" %(count))
                break
    return ClusterList

def IdentifyVgeneCluster(sMat):
    ## Input Vgene score matrix
    vG={}
    n=len(sMat)
    IDs=[x for x in range(n)]
    for kk in IDs:
        LL=sMat[:,kk]
        vL=np.where(LL>=thr_v)[0]
        if len(vL)>0:
            vG[kk]=vL
    CL=IdentifyMotifCluster(vG)
    return CL
    
def ParseFa(fname):
    InputStr=open(fname).readlines()
    FaDict={}
    seq=''
    for line in InputStr:
        if line.startswith('>'):
            if len(seq)>0:
                FaDict[seqHead]=seq
                seq=''
            seqHead=line.strip()
        else:
            seq+=line.strip()
    if seqHead not in FaDict:
        FaDict[seqHead]=seq
    return FaDict

def PreCalculateVgeneDist(VgeneFa="Imgt_Human_TRBV.fasta"):
    ## Only run one time if needed
    FaDict=ParseFa(cur_dir+VgeneFa)
    VScore={}
    CDR1Dict={}
    CDR2Dict={}
    for kk in FaDict:
        if '|' in kk:
            VV=kk.split('|')[1]
        else:
            VV=kk[1:]
        CDR1Dict[VV]=FaDict[kk][26:37]  ## Imgt CDR1: 27 - 38
        CDR2Dict[VV]=FaDict[kk][55:64]  ## Imgt CDR2: 56 - 65
    Vkeys=list(CDR1Dict.keys())
    nn=len(Vkeys)
    for ii in range(0,nn):
        V1=Vkeys[ii]
        s1_CDR1=CDR1Dict[V1]
        s1_CDR2=CDR2Dict[V1]
        for jj in range(ii,nn):
            V2=Vkeys[jj]
            s2_CDR1=CDR1Dict[V2]
            s2_CDR2=CDR2Dict[V2]
            score1=SeqComparison(s1_CDR1,s2_CDR1)
            score2=SeqComparison(s2_CDR2,s2_CDR2)
            #print score1+score2
            VScore[(V1,V2)]=score1+score2
    gg=open('VgeneScores.txt','w')
    for kk in VScore:
        vv=VScore[kk]
        line=kk[0]+'\t'+kk[1]+'\t'+str(vv)+'\n'
        gg.write(line)
    gg.close()

def EncodeRepertoire(inputfile, outdir, outfile='',exact=True, ST=3, thr_v=3.7, thr_s=3.5, VDict={},Vgene=True,thr_iso=10, gap=-6, GPU=False,Mat=False, verbose=False):
    ## No V gene version
    ## Encode CDR3 sequences into 96 dimensional space and perform k-means clustering
    ## If exact is True, SW alignment will be performed within each cluster after isometric encoding and clustering
    h=open(inputfile)
    t1=time.time()
    alines=h.readlines()
    ww=alines[0].strip().split('\t')
    if not ww[0].startswith('C'):
        ## header line
        hline=alines[0]
        alines=alines[1:]        
    elif 'CDR3' in ww[0]:
        hline=alines[0]
        alines=alines[1:]
    else:
        hline='CDR3\t'+'\t'.join(['Info'+str(x) for x in range(len(ww)-1)])        
    seqs=[]
    vgs=[]
    infoList=[]
    count=0
    if verbose:
        print('Creating CDR3 list')
    for ll in alines:
        ww=ll.strip().split('\t')
        cdr3=ww[0]
        if '*' in cdr3:
            continue
        if '_' in cdr3:
            continue
        seqs.append(ww[0])
        if Vgene:
            vgs.append(ww[1])
            infoList.append('\t'.join(ww[1:]))
        else:
            infoList.append('\t'.join(ww[1:]))
        count+=1
    if len(outfile)==0:
        outfile=inputfile.split('/')
        outfile=outfile[len(outfile)-1]
        outfile=outdir+'/'+re.sub('\\.[txcsv]+','',outfile)+'-'+'-RotationEncodingBL62.txt'
    g=open(outfile,'w')
    tm=strftime("%Y-%m-%d %H:%M:%S", gmtime())
    InfoLine='##TIME:'+tm+'|cmd: '+sys.argv[0]+'|'+inputfile+'|IsometricDistance_Thr='+str(thr_iso)+'|thr_v='+str(thr_v)+'|thr_s='+str(thr_s)+'|exact='+str(exact)+'|Vgene='+str(Vgene)+'|ST='+str(ST)
    g.write(InfoLine+'\n')
    g.write("##Column Info: CDR3 aa sequence, cluster id, other information in the input file\n")
    gr=0
    ## Split into different lengths
    LD,VD, ID,SD= BuildLengthDict(seqs, vGene=vgs,INFO=infoList,sIDs=[x for x in range(len(seqs))])
    LDu, VDu, IDu, SDu = CollapseUnique(LD, VD, ID, SD)
    if Mat:
        Mfile=outfile+'_EncodingMatrix.txt'
        h=open(Mfile, 'w')
    for kk in LDu:
        if verbose:
            print("---Process CDR3s with length %d ---" %(kk))
        vSD=LDu[kk]
        vSD0=[x for x in range(len(vSD))]
        vss=SDu[kk]
        vInfo=IDu[kk]
        flagL=[len(x)-1 for x in vInfo]
        if verbose:
            print(' Performing CDR3 encoding')
        dM=np.array([EncodingCDR3(x[ST:-2], M6, n0) for x in vss])
        dM=dM.astype("float32")
        if verbose:
            print(" The number of sequences is %d" %(dM.shape[0]))
        if Mat:
            for ii in range(len(vss)):
                line=vss[ii]+'\t'+vInfo[ii][0]+'\t'
                NUMs=[str(xx) for xx in dM[ii,:]]
                line += '\t'.join(NUMs) + '\n'
                h.write(line)
        sID=[x for x in range(dM.shape[0])]
        t2=time.time()
        if verbose:
            print(' Done! Total time elapsed %f' %(t2-t1))
        Cls = ClusterCDR3(dM, flagL, thr=thr_iso - 0.5*(15-kk), verbose=verbose)  ## change cutoff with different lengths
        if verbose:
            print("     Handling identical CDR3 groups")
        Cls_u=[]
        for ii in range(len(Cls)):
            cc=Cls[ii]
            if len(cc) == 1:
                ## Handle identical CDR3 groups first
                if flagL[cc[0]]>0:
                    gr += 1
                    jj=cc[0]
                    for v_info in vInfo[jj]:
                        line=vss[jj]+'\t'+str(gr)+'\t'+v_info+'\n'
                        _=g.write(line)
            else:
                Cls_u.append(cc)
        Cls=Cls_u
        t2=time.time()
        if verbose:
            print(' Done! Total time elapsed %f' %(t2-t1))
        if Vgene:
            vVgene=VDu[kk]
            if verbose:
                print('     Matching variable genes')
            Cls_v=[]
            for cc in Cls:
                Nc=len(cc)
                sMat={}
                for ii in range(Nc):
                    v1=vVgene[cc[ii]]
                    for jj in range(ii,Nc):
                        if jj==ii:
                            continue
                        v2=vVgene[cc[jj]]
                        if (v1, v2) not in VDict:
                            if v1 == v2:
                                if ii not in sMat:
                                    sMat[ii]=[jj]
                                else:
                                    sMat[ii].append(jj)
                                if jj not in sMat:
                                    sMat[jj]=[ii]
                                else:
                                    sMat[jj].append(ii)
                            continue
                        if VDict[(v1,v2)] >= thr_v:
                                if ii not in sMat:
                                    sMat[ii]=[jj]
                                else:
                                    sMat[ii].append(jj)
                                if jj not in sMat:
                                    sMat[jj]=[ii]
                                else:
                                    sMat[jj].append(ii)
                vCL=IdentifyMotifCluster(sMat)
                vCL_List=list(chain(*vCL))
                for ii in range(Nc):
                    uu=flagL[cc[ii]]
                    if uu>0 and ii not in vCL_List:
                        vCL.append([ii])
                for vcc in vCL:
                    Cls_v.append(list(np.array(cc)[np.array(vcc)]))
            Cls=[]
            for ii in range(len(Cls_v)):
                cc=Cls_v[ii]
                if len(cc) == 1:
                    ## Handle identical CDR3 groups first
                    gr += 1
                    jj=cc[0]
                    for v_info in vInfo[jj]:
                        line=vss[jj]+'\t'+str(gr)+'\t'+v_info+'\n'
                        _=g.write(line)
                else:
                    Cls.append(cc)
        if exact:
            if verbose:
                print(' Performing Smith-Waterman alignment')
            Cls_s=[]
            for cc in Cls:
                Nc=len(cc)
                if len(cc)<=3:
                    sMat=np.zeros((Nc,Nc))
                    for ii in range(Nc):
                        s1=vss[cc[ii]]
                        for jj in range(ii,Nc):
                            if jj==ii:
                                continue
                            s2=vss[cc[jj]]
                            if len(s1) != len(s2):
                                continue
                            if len(s1)<=5:
                                continue
                            sw=SeqComparison(s1[ST:-2],s2[ST:-2],gap=gap)
                            sw=sw/(len(s1)-ST-2)
                            sMat[ii,jj]=sw
                            sMat[jj,ii]=sw
                    s_max=[]
                    for ii in range(Nc):
                        s_max.append(np.max(sMat[:,ii]))
                    cc_new=[]
                    for ii in range(Nc):
                        if s_max[ii]>=thr_s:
                            cc_new.append(cc[ii])
                    if len(cc_new)>1:
                        Cls_s.append(cc_new)
                    else:
                        for ii in range(Nc):
                            uu=flagL[cc[ii]]
                            if uu>0:
                                Cls_s.append([cc[ii]])
#                    print(Cls_s)
                    Cls_sList=list(chain(*Cls_s))
                    for ii in range(len(cc)):
                        uu=flagL[cc[ii]]
                        if uu>0 and cc[ii] not in Cls_sList:
                            Cls_s.append([cc[ii]])
                else:
                    CDR3s=[vss[x] for x in cc]
                    sIDs=np.array([vSD0[x] for x in cc])
                    sIDs0=[x for x in range(len(cc))]
                    Kset=KmerSet(CDR3s, sIDs0, KS=5, st=ST, ed=2)
                    SSG=generateSSG(Kset, CDR3s, k_thr=1)
                    tmpVgenes=['TRBV2']*len(CDR3s)
                    SSGnew=UpdateSSG(SSG, CDR3s, tmpVgenes, Vscore=VDict, cutoff=thr_s+4)
                    CLall=IdentifyMotifCluster(SSGnew)
                    CLall_list=list(chain(*CLall))
                    for ii in range(len(cc)):
                        uu=flagL[cc[ii]]
                        if uu>0 and ii not in CLall_list:
                            CLall.append([ii])
                    for cl in CLall:
                        ccs=list(sIDs[np.array(cl)])
                        Cls_s.append(ccs)
            Cls=Cls_s
        if verbose:
            print(' Writing results into file')
        for ii in range(len(Cls)):
#            if ii % 100000 == 0 and ii>0:
                #print('      %d sequences written' %(ii))
            cc=Cls[ii]
            gr+=1
            for jj in cc:
                for v_info in vInfo[jj]:
                    line=vss[jj]+'\t'+str(gr)+'\t'+v_info+'\n'
                    _=g.write(line)
    g.close()
    if Mat:
        h.close()

def OrderUnique(Ig):
    vv=list(Ig.values())
    kk=list(Ig.keys())
    LL=[len(x[1]) for x in vv]
    v0=[x[0][0] for x in vv]
    v1=[x[0][1] for x in vv]
    zkk=zip(kk,v0,v1,LL)
    zkks=sorted(zkk,key=lambda x: (x[1],x[3]))
    nk=len(zkks)
    keep_id=[0]
    ii=1
    n_pre=str(zkks[0][1])+'_'+str(zkks[0][2])
    while ii<nk:
        n_cur=str(zkks[ii][1])+'_'+str(zkks[ii][2])
        if n_cur==n_pre:
            ii+=1
            continue
        else:
            keep_id.append(ii)
            n_pre=n_cur
            ii+=1
            continue
    nid=[x[0] for x in zkks]
    filtered_id=np.array(nid)[np.array(keep_id)]
    Igs={}
    for ii in filtered_id:
        Igs[kk[ii]]=vv[ii]
    return Igs, filtered_id

def ClusterCDR3(dM, flagL, thr=10, GPU=False, verbose=False):
    ## flagL: flag vector for identical CDR3 groups, >0 for grouped non-identical CDR3s
    Cls=[]
    flag=0
    dM1=dM
    flagL=np.array(flagL)
    if GPU:
        res = faiss.StandardGpuResources()
    while 1:
#        print("     %d number of clusters, with %d sequences" %(len(Cls),dM1.shape[0]))
        if verbose:
            print('=',end='')
        index = faiss.IndexFlatL2(Ndim*6)
        if GPU:
            index = faiss.index_cpu_to_gpu(res, 0, index)
        index.add(dM1)
        if flag==0:
            D, I = index.search(dM1, 2)
            vv=np.where((D[:,1]<=thr))[0]
            vv0=np.where((D[:,1]>thr) & (flagL>0))[0]
            for v in vv0:
                Cls.append([v])
            tmp_dM=np.zeros((len(vv),Ndim*6))
            Ig_new={}
            for ii in range(len(vv)):
                v=vv[ii]
                Idx=I[v,]
                if v not in Idx:
                    Idx[0]=v
                Ig_new[ii]=(sorted(list(set(Idx))),sorted(list(set(Idx))))
                tmp_dM[ii,]=(dM1[Idx[0],]+dM1[Idx[1],])/2
            if len(Ig_new)==0:
                if verbose:
                    print('type 0 break')
                break
#                print('%d of sequence left at cycle %d' %(len(Ig_new),flag))
            Igs, fid=OrderUnique(Ig_new)
            tmp_dM=tmp_dM[fid,]
            Ig_new=Igs
        else:
            D, I = index.search(dM1,2)
            vv=np.where(D[:,1]<=thr)[0]
            vv0=np.where(D[:,1]>thr)[0]
            ## move groups in vv0 to Cls
            kkg=list(Ig.keys())
            for v in vv0:
                ng=list(Ig[kkg[v]][1])
    #            if ng not in Cls:
                Cls.append(ng)
            tmp_dM=np.zeros((len(vv),Ndim*6))
            Ig_new={}
            for ii in range(len(vv)):
                v=vv[ii]
                idx1=I[v,0]
                idx2=I[v,1]
                if v not in I[v,]:
                    idx1=v
#                Ig_new[ii]=sorted(list(set(list(Ig[kkg[idx1]])+list(Ig[kkg[idx2]]))))
                Ig_new[ii]=(sorted(list(set([idx1,idx2]))),  ## First entry records the relative index of a sequence clique
                            sorted(list(set(list(Ig[kkg[idx1]][1])+list(Ig[kkg[idx2]][1])))))  ## Second entry records the absolute index of a sequence
                tmp_dM[ii,]=(dM1[idx1,]+dM1[idx2,])/2
            if len(Ig_new)==0:
                if verbose:
                    print("\ntype I break")
                kkg=list(Ig.keys())
                for kk in kkg:
                    ng=list(Ig[kk][1])
                    if ng not in Cls:
                        Cls.append(ng)
                break
#            print('%d of sequence left at cycle %d' %(len(Ig_new),flag))
            Igs, fid=OrderUnique(Ig_new)
            tmp_dM=tmp_dM[fid,]
            Ig_new=Igs
        if flag>0:
            if Ig == Ig_new:
                if verbose:
                    print("\ntype II break")
                kkg=list(Ig.keys())
                for kk in kkg:
                    ng=list(Ig[kk][1])
                    if ng in Cls:
                        continue
                    Cls.append(ng)
                break
        Ig=Ig_new
        tmp_dM=tmp_dM.astype('float32')
        dM1=tmp_dM
        flag+=1
    return Cls

def ClusterCDR3r(dM, flagL, thr = 10, verbose = False):
    index = faiss.IndexFlatL2(Ndim*6)
    index.add(dM)
    lims, D, I = index.range_search(dM, thr)
    # with open('cdr3.npy', 'wb') as f:
    #     np.save(f, lims)
    #     np.save(f, D)
    #     np.save(f, I)
    #     np.save(f, dM)
    
    # now clustering results
    N = dM.shape[0]
    neighborSize = np.array([lims[cur_idx_i+1] - lims[cur_idx_i] for cur_idx_i in range(N)])
    # to_cluster = np.ones( (N,))
    clusterNo = 0
    cluster = - np.ones( (N, ),  dtype = np.int32)
    idx = np.where(cluster < 0)[0]
    unclustered = [np.argmax(neighborSize[idx])]
    depth = 0
    while True:
        if len(unclustered) == 0: break
        # cur_idx = unclustered[0] # first unclustered index
        cur_idx = unclustered
        cluster[cur_idx] = clusterNo # assign cluster
   
        neighbor = np.unique(np.array(list(chain (* [I[(lims[cur_idx_i]): lims[cur_idx_i+1]] for cur_idx_i in cur_idx]))))
        # find those unclusterred
        idx = np.where(cluster[neighbor] < 0)[0]
        if len(idx) == 0:
            depth = 0
            clusterNo += 1
            idx = np.where(cluster < 0)[0]
            if len(idx) == 0: break
            unclustered = [idx[np.argmax(neighborSize[idx])]]
           
        else:
            if depth > 3:
                depth = 0
                clusterNo += 1
            unclustered = neighbor[idx]
            depth += 1
#    print('clusterNo = ', clusterNo)
    Cls = [ [] for i in range(clusterNo)]
    for idx, i in enumerate(cluster):
            Cls[i].append(idx)
#    print("Cls[:5] = ", Cls[:5])
#    print("len(Cls) = ", len(Cls),
#          ', #elem=', sum([len(i) for i in Cls]),
#          ', #single=', sum([len(i) for i in Cls if len(i) == 1]),
#          ', #non_single=', sum([len(i) for i in Cls if len(i) != 1]),          
#          ', #max=', max([len(i) for i in Cls]))
    return Cls

def CommandLineParser():
    parser=OptionParser()
    print ('''
GIANA: Geometric Isometry based ANtigen-specific tcr Alignment
Ultrafast short peptide alignment exclusively designed for large-scale adaptome analysis

Input columns:
1. CDR3 amino acid sequence (Starting from C, ending with the first F/L in motif [FL]G.G)
2. Variable gene name in Imgt format: TRBVXX-XX*XX
3. Joining gene name (optional)
4. Frequency (optional)
5. Other information (optional)

!!! ALL amino acid letters must be CAPITAL !!!

''')
    parser.add_option("-d","--directory",dest="Directory",help="Input repertoire sequencing file directory. Please make sure that all the files in the directory are input files.",default="")
    parser.add_option("-f","--file",dest="File",default='',help="Input single file of CDR3 sequences for grouping")
    parser.add_option("-F","--fileList",dest="files",default='',help='Alternative input: a file containing the full path to all the files. If given, overwrite -d and -f option')
    parser.add_option("-t","--threshold",dest="thr",default=7,help="Isometric distance threshold for calling similar CDR3 groups. Without -E, smaller value will increase speed. With -E, smaller value will increase specificity. Must be smaller than 12.")
    parser.add_option("-S","--threshold_score",dest="thr_s",default=3.6, help="Threshold for Smith-Waterman alignment score (normalized by CDR3 length). Default 3.6")
    parser.add_option("-G","--threshold_vgene",dest="thr_v",default=3.7,help="Threshold for variable gene comparison. Default 3.7.")
    parser.add_option("-o","--output",dest="OutDir",default='./',help="Output directory for intermediate and final outputs.")
    parser.add_option("-O","--outfile",dest="OutFile",default='',help="Output file name. If not given, a file with --RotationEncoding will be added to the input file as the output file name.")
    parser.add_option("-T","--startPosition",dest='ST',default=3, help="Starting position of CDR3 sequence. The first ST letters are omitted. CDR3 sequence length L must be >= ST+7 ")
    parser.add_option("-g","--GapPenalty",dest="Gap",default= -6,help="Gap penalty,default= -6. Not used.")
    parser.add_option("-n","--GapNumber",dest="GapN",default=1,help="Maximum number of gaps allowed when performing alignment. Max=1, default=1. Not used.")
    parser.add_option("-V","--VariableGeneFa",dest="VFa",default="Imgt_Human_TRBV.fasta",help="IMGT Human beta variable gene sequences")
    parser.add_option("-v","--VariableGene",dest="V",default=True,action="store_false",help="If False, GIANA will omit variable gene information and use CDR3 sequences only. This will yield reduced specificity. The cut-off will automatically become the current value-4.0")
    parser.add_option("-e","--Exact",dest="E",default=True,action="store_false",help="If False, GIANA will not perform Smith-Waterman alignment after isometric encoding.")
    parser.add_option("-N","--NumberOfThreads",dest="NN",default=1,help="Number of threads for multiple processing. Not working so well.")
    parser.add_option("-M","--EncodingMatrix", dest="Mat", default=False,action="store_true", help="If true, GIANA will export the isometric encoding matrix for each TCR. Default: False.")
    parser.add_option("-U","--UseGPU",dest="GPU", default=False, action="store_true",help="Use GPU for Faiss indexing. Must be CUDA GPUs.")
    parser.add_option("-q","--queryFile",dest="Query",default='',help="Input query file, if given, GIANA will run in query mode, also need to provide -r option.")
    parser.add_option("-r","--refFile",dest="ref", default='',help="Input reference file. Query model required.")
    parser.add_option("-b","--Verbose", dest='v', default=False, action="store_true", help="Verbose option: if given, GIANA will print intermediate messages.")
    return parser.parse_args()

def main():
    (opt,_)=CommandLineParser()
    cutoff=float(opt.thr)
    OutDir=opt.OutDir
    thr_s=float(opt.thr_s)
    ## Check if query mode first
    qFile=opt.Query
    if len(qFile)>0:
        ## query mode
        t1=time.time()
        if qFile.endswith('/'):
            ## input query is a directory
            qFs=os.listdir(qFile)
            qFileList=[]
            for ff in qFs:
                qFileList.append(qFile+ff)
        else:
            qFileList=[qFile]
        rFile=opt.ref
        if len(rFile)==0:
            raise("Must provide reference file in query mode!")
        else:
            ## check if reference cluster file exists
            rFile0=re.sub('\\.txt','',rFile)
            refClusterFile=rFile0+'--RotationEncodingBL62.txt'
            if not os.path.exists(refClusterFile):
                raise("Must run clustering on reference file first! Did you forget to put the clustering file in this directory?")
            rData=CreateReference(rFile)
            t2=time.time()
            print("Reference created. Elapsed %f" %(t2-t1))
            for qf in qFileList:
                t2_0=time.time()
                print("Querying "+qf)
                qf_s=qf.split('/')[-1]
                #outFile=re.sub('\\.txt','',qf_s)+'_query_'+rFile0+'.txt'
                outFile=os.path.splitext(qf_s)[0]+'_query_'+os.path.basename(rFile0)+'.txt'
                of=OutDir+'/'+outFile
                if path.exists(of):
                    print(of+' already exits. Skipping.')
                    continue
                MakeQuery(qf, rData, thr=cutoff, thr_s=thr_s)
                t2=time.time()
                print("     Build query clustering file. Elapsed %f" %(t2-t1))
                print("Now mering with reference cluster")
                MergeExist(refClusterFile, OutDir+'/'+outFile)
                t2=time.time()
                print("         Time of elapsed for query %s: %f" %(qf, t2-t2_0))
    else:
        ## regular clustering mode
        FileDir=opt.Directory
        if len(FileDir)>0:
                files=os.listdir(FileDir)
                files0=[]
                for ff in files:
                        ff=FileDir+'/'+ff
                        files0.append(ff)
                files=files0
        else:
                files=[]
        File=opt.File
        if len(File)>0:
                files=[File]
        FileList=opt.files
        if len(FileList)>0:
                files=[]
                fL=open(FileList)
                for ff in fL.readlines():
                        files.append(ff.strip())
        VFa=opt.VFa
        PreCalculateVgeneDist(VFa)
        vf=open('./VgeneScores.txt')  ## Use tcrDist's Vgene 80-score calculation
        VScore={}
        VV=opt.V
        EE=opt.E
        Mat=opt.Mat
        ST=int(opt.ST)
        thr_v=float(opt.thr_v)
        verbose=opt.v
        if VV:
            while 1:
                line=vf.readline()
                if len(line)==0:
                    break
                ww=line.strip().split('\t')
                VScore[(ww[0],ww[1])]=int(float(ww[2]))/20
                VScore[(ww[1],ww[0])]=int(float(ww[2]))/20
        Gap=int(opt.Gap)
        Gapn=int(opt.GapN)
        OutFile=opt.OutFile
        GPU=opt.GPU
        st=3
        ed=1
        NT=int(opt.NN)
        faiss.omp_set_num_threads(NT)
        for ff in files:
            print("Processing %s" %ff)
            EncodeRepertoire(ff, OutDir, OutFile, ST=ST, thr_s=thr_s, thr_v=thr_v, exact=EE,VDict=VScore, Vgene=VV, thr_iso=cutoff, gap=Gap, GPU=GPU, Mat=Mat, verbose=verbose)
        
if __name__ == "__main__":
    t0=time.time()
    main()
    print ("Total time elapsed: %f" %(time.time()-t0))
    print ("Maximum memory usage: %f MB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000))

