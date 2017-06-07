

function   R= Fun_PCADimReduction (FV , MaskAct , InfThr)

[FV_H,FV_W , FV_num ] =size((FV)) ;
FV=double(FV) ;


% get MaskActivity ----------------------
MaskActLevel= imresize(MaskAct,[FV_H,FV_W]);
MaskActLevel(MaskActLevel>=4.9)=5 ; MaskActLevel(MaskActLevel<=1.1)=1 ;
MaskActLevel((MaskActLevel<5)&(MaskActLevel>1))=3 ;
MaskALevel4F =(MaskActLevel==1) ;
MaskALevel4B =(MaskActLevel==5) ;
MaskALevel4U =(MaskActLevel==0) ;
MaskALevel4A =(MaskActLevel==3) ;

%----------------------------------------


FV_Vect = reshape(FV,[FV_H*FV_W,FV_num]);
FV=[] ;

%S1) Get PCA of Data ==================================================
MaskALevelTrain = (MaskALevel4F | MaskALevel4B) ;
TCDataPCATrainB=[] ; TCDataPCATrainF=[] ; TCDataPCATrain=[] ;
for i=1 :FV_num
    
    a= FV_Vect(:,i);
    
    
    TCDataPCATrainF (:,i) = a(MaskALevel4F(:));
    TCDataPCATrainB (:,i) = a(MaskALevel4B(:));
end
a=[] ;
%% push F and B sets to have equal number of samples ---------------

[FSampleNum,DimNumB4PCA ] = size(TCDataPCATrainF);
BSampleNum = size(TCDataPCATrainB,1);
TCDataPCATrainF_Ext= TCDataPCATrainF ;
TCDataPCATrainB_Ext= TCDataPCATrainB ;
if FSampleNum >=  BSampleNum % add sample to backgrounds
    
    NumAddSample = round(abs(FSampleNum-BSampleNum)) ;
    TInd= round(rand(NumAddSample,1)*(BSampleNum-1))+1 ;
    TSample = TCDataPCATrainB(TInd,:)+rand(NumAddSample,DimNumB4PCA) ; % new samples
    TCDataPCATrainB_Ext(BSampleNum+1 :BSampleNum+NumAddSample,:)=TSample ; % extended set
else  % add samples to foreground set
    NumAddSample = round(abs(BSampleNum-FSampleNum)) ;
    TInd= round(rand(NumAddSample,1)*(FSampleNum-1))+1 ;
    TSample = TCDataPCATrainF(TInd,:)+rand(NumAddSample,DimNumB4PCA) ; % new samples
    TCDataPCATrainF_Ext(FSampleNum+1 :FSampleNum+NumAddSample,:)=TSample ; % extended set
end % eof if
TSample=[] ; TInd=[] ;
TCDataPCATrainF=[] ; TCDataPCATrainB=[] ;

TCTrainSet_Ext = [TCDataPCATrainF_Ext ; TCDataPCATrainB_Ext]; % Train Extended Data

TCDataPCATrainF_Ext= [] ; TCDataPCATrainF_Ext=[] ;



PCA_InoThreshold = InfThr ;  % information threshold to select best eigen vectors

V=[] ;
%         V= cov(TCData) ;
V= cov(TCTrainSet_Ext) ; % use extended sets
[TPCA_Coef,TPCA_Latent,TPCA_Score] = pcacov((V)) ;
sel=1 ; % find top fist eign vector whic keep .9 percent of informations
while sum(TPCA_Latent(1:sel))/sum(TPCA_Latent)  <=PCA_InoThreshold
    sel=sel+1  ;
end

TCDataPCA =   FV_Vect*TPCA_Coef(:,1:sel) ; % FV in PCA Format -------
TCDataPCA4F=[] ; TCDataPCA4B=[] ; TCDataPCA4U=[] ; TCDataPCA4A=[] ;



%% Return Result --------------------------------------------------------


R = reshape(TCDataPCA,[FV_H,FV_W,sel]) ; % reduced dimension of data ;






end
