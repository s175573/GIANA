## Query mode of GIANA, to be loaded by GIANA, cannot run alone

import shelve
import subprocess as sp
import pandas as pd
from GIANA4 import *

def CreateReference(rFile, outdir='./', Vgene=True, ST=3):
    ## convert input reference file into a python workplace
    h=open(rFile)
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
    LD,VD, ID,SD= BuildLengthDict(seqs, vGene=vgs,INFO=infoList,sIDs=[x for x in range(len(seqs))])
    LDu_r, VDu_r, IDu_r, SDu_r = CollapseUnique(LD, VD, ID, SD)
    flagLD_r={}
    dMD_r={}
    for kk in LDu_r:
        vss=SDu_r[kk]
        vInfo=IDu_r[kk]
        flagL=[len(x)-1 for x in vInfo]
        flagLD_r[kk]=flagL
        dM=np.array([EncodingCDR3(x[ST:-2], M6, n0) for x in vss])
        dM=dM.astype("float32")        
        dMD_r[kk]=dM
##    ff0=re.sub('.txt','',rFile)
##    outfile=outdir+ff0+'_giana_ref.shelve'
##    giana_shelf = shelve.open(outfile, 'n')
##    giana_shelf['flagLD']=flagLD_r
##    giana_shelf['dMD']=dMD_r
##    giana_shelf['LDu']=LDu_r
##    giana_shelf['VDu']=VDu_r
##    giana_shelf['IDu']=IDu_r
##    giana_shelf['SDu']=SDu_r
##    giana_shelf.close()
    return [LDu_r, VDu_r, IDu_r, SDu_r, dMD_r]

def MakeQuery(qFile, rData=[],dbFile=None, Vgene=True, thr=7, ST=3, thr_s=3.3):
    if dbFile is not None:
        with shelve.open(dbFile) as db:
            for key in db:
                globals()[key]=db[key]
    else:
        if len(rData)==0:
            raise("Need to provide either a reference file or a shelve")
        LDu_r=rData[0]
        VDu_r=rData[1]
        IDu_r=rData[2]
        SDu_r=rData[3]
        dMD_r=rData[4]
    h=open(qFile)
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
    LD,VD, ID,SD= BuildLengthDict(seqs, vGene=vgs,INFO=infoList,sIDs=[x for x in range(len(seqs))])
    LDu, VDu, IDu, SDu = CollapseUnique(LD, VD, ID, SD)
    tmpFile='tmp_query.txt'
    g=open(tmpFile,'w')
    for kk in LDu:
        vss=SDu[kk]
        vInfo=IDu[kk]
        vss_r=SDu_r[kk]
        vInfo_r=IDu_r[kk]
        flagL=[len(x)-1 for x in vInfo]
        dM_r=dMD_r[kk]
        dM=np.array([EncodingCDR3(x[ST:-2], M6, n0) for x in vss])
        dM=dM.astype("float32")
        nq=dM.shape[0]
        nr=dM_r.shape[0]
        vssc=vss+vss_r
        vInfoc=vInfo+vInfo_r
        dMc=np.concatenate((dM, dM_r))
        index = faiss.IndexFlatL2(Ndim*6)
        index.add(dMc)
        D, I = index.search(dM, 2)
        vv=np.where((D[0:nq,1]<=thr))[0]
        flagL=np.array(flagL)
        vv0=np.where((D[0:nq,1]>thr) & (flagL>0))[0]
        curList=[]
        for v in vv0:
            for ii in range(len(vInfoc[v])):
                line=vssc[v]+'\t'+vInfoc[v][ii]+'\t'+'query\n'
                _=g.write(line)
        for v in vv:
            tmpI=I[v,]
            if v not in tmpI:
                tmpI[0]=v
            idx1=tmpI[0]
            idx2=tmpI[1]
            c1=vssc[idx1]
            c2=vssc[idx2]
            info1=vInfoc[idx1]
            info2=vInfoc[idx2]
            for tmpInfo in info1:
                tup1=(c1, tmpInfo)
                if tup1 not in curList:
                    if idx1<nq:
                        line1=c1+'\t'+tmpInfo+'\t'+'query\n'
                    else:
                        line1=c1+'\t'+tmpInfo+'\t'+'ref\n'
                    _=g.write(line1)
                    curList.append(tup1)
            for tmpInfo in info2:
                tup2=(c2, tmpInfo)
                if tup2 not in curList:
                    if idx2<nq:
                        line2=c2+'\t'+tmpInfo+'\t'+'query\n'
                    else:
                        line2=c2+'\t'+tmpInfo+'\t'+'ref\n'
                    _=g.write(line2)
                    curList.append(tup2)
    g.close()
    cmd='python3 GIANA4.py -f tmp_query.txt -S '+str(thr_s)
    p=sp.run(cmd, shell=True)

def MergeExist(refClusterFile, outFile='queryFinal.txt',queryClusterFile='tmp_query--RotationEncodingBL62.txt', direction='q'):
    ## This function compare the query file with ref cluster file and merge the two based on shared TCRs
    ## If direction is 'q', the overlapping clusters will be added to the query file
    ## If direction is 'r', the overlapping and non-overlapping clusters will be added to the reference file
    refT=pd.read_table(refClusterFile, skiprows=2, delimiter='\t', header=None)
    queryT=pd.read_table(queryClusterFile, skiprows=2, delimiter='\t', header=None)
    nq=queryT.shape[1]
    nr=refT.shape[1]
    if nr != nq-1:
        print("ERROR: Make sure reference and the query samples have the same columns!")
        print("No query file is generated.")
        return
    gn=np.unique(queryT[1])
    queryTs=pd.DataFrame([], columns=queryT.columns)
    for nn in gn:
            tmp_ddq=queryT.loc[np.where(queryT[1]==nn)[0],:]
            cls_lab=np.unique(tmp_ddq[nq-1])
            if len(cls_lab)==1:
                if cls_lab[0]=='ref':
                    continue
            queryTs=queryTs._append(tmp_ddq)
    queryTs.index=range(queryTs.shape[0])
    keyr=refT[0]+'_'+refT[2]
    keyq=queryTs[0]+'_'+queryTs[2]
    vvr=np.where(queryTs[nq-1]=='ref')[0]
    vvr_in=np.where(keyr.isin(keyq[vvr]))[0]
    gn_r=list(refT.loc[vvr_in,1].drop_duplicates())
    ddo=pd.DataFrame([], columns=refT.columns)
    for nn in gn_r:
        tmp_dd=refT.loc[np.where(refT[1]==nn)[0],:]
        tmpkey=tmp_dd[0]+'_'+tmp_dd[2]
        vv=np.where(keyq.isin(tmpkey))[0][0]
        gq=queryTs[1][vv]
        tmp_dd[1]=gq
        ddo=ddo._append(tmp_dd)
    if direction=='q':
        ddo[nq-1]='ref'
        ## remove groups that contain only ref group
        queryTs=queryTs._append(ddo)
        queryTs=queryTs.drop_duplicates()
        queryTs.to_csv(outFile, sep='\t',header=False,index=False)
#    queryTs.index=range(queryTs.shape[0])
    if direction=='r':
        ## to be developed
        pass












