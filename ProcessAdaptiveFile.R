ProcessAdaptiveVgenes <- function(dd){
  Vgene_o=c(1,2,9,13:19,26:28,30)
  Vgene_o=as.character(Vgene_o)
  dd=dd[which(nchar(as.character(dd[,2]))>0),]
  gsub('TCRBV[0]{0,1}','TRBV',dd[,2])->tmpV
  ## Multiple calls
  vv.m=grep('/',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],'/')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  vv.m=grep(',',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],',')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  tmpV=gsub('-0','-',tmpV)
  v_digit=grep('\\*',tmpV)
  if(length(v_digit)>0)tmpV[- v_digit] = paste(tmpV[- v_digit],'*01',sep='') else tmpV=paste(tmpV,'*01',sep='')
  ## 1. Orphan V genes do not have "-1", need to remove
  Vnumbers=gsub('TRBV','',tmpV)
  Vnumbers=gsub('\\*.+','',Vnumbers)
  Vnumbers=gsub('-.+','',Vnumbers)
  vv.o=which(Vnumbers %in% Vgene_o)
  tmpV1=tmpV
  tmpV1[vv.o]=gsub('-1','',tmpV1[vv.o])
  ## 2. Non-orphan V genes but without "-1", need to add
  vv.no=which(! Vnumbers %in% Vgene_o)
  vv.non=grep('-',tmpV1[vv.no])
  if(length(vv.non)>0)tmpV1[vv.no][-vv.non]=gsub('\\*01','-1*01',tmpV1[vv.no][-vv.non])
  dd[,2]=tmpV1
  return(dd)
}


PrepareAdaptiveFile <- function(indir,outdir,thr=10000){
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('\\.tsv|\\.txt',ff))==0 ) next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    print(ff)
    ddnew=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    tmp.vv.rm=grep('\\*',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.vv.rm=grep('X',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.nn=nchar(ddnew[,2])
    tmp.vv=which(tmp.nn>=10 & tmp.nn<=24)
    #tmp.vv=which(tmp.nn>=1)
    ddnew=ddnew[tmp.vv,]
    ddnew=ddnew[which(ddnew[,6]!='unresolved'),]
    ddnew=ddnew[grep('^C.+[FV]$',ddnew[,2]),]
    #ss=sum(ddnew[,3])
    #ddnew[,4]=ddnew[,3]/ss
    ddnew=ddnew[order(ddnew[,4],decreasing=T),]
    #q75=quantile(ddnew[,3],0.75)
    #THR=min(ddnew[,3]+1)
    #tmp.vv=which(ddnew[,3]>max(q75,THR))
    #if(length(tmp.vv)<=thr)tmp.vv=1:thr
    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c(2,6,4)] else ddnew=ddnew[,c(2,6,4)]
    #ddnew=ddnew[1:5000,]
    ddnew=ddnew[1:min(thr,nrow(ddnew)),]
    ddnew=ProcessAdaptiveVgenes(ddnew)
    ddnew=cbind(ddnew, RANK=rank(ddnew[,3])/nrow(ddnew))
    TMP=gsub('.+/TCR_peptide/data/Adaptive/','',indir)
    tmp=strsplit(TMP,'/')[[1]][[1]]
    ddnew=cbind(ddnew, Info=paste(tmp,ff0,sep=':'))
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }
}

