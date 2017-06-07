
function   [RFVLDA , RLDA_Eig ]= Fun_LDADimReduction (FV , MaskAct , InfThr,FunMod )



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

for i=1 : FV_num
    a= FV_Vect(:,i);
    TCData4F (:,i) = a(MaskALevel4F(:));
    TCData4B (:,i)= a(MaskALevel4B(:));
    TCData4U (:,i)= a(MaskALevel4U(:));
    TCData4A(:,i) = a(MaskALevel4A(:));
end
a= [] ;

%% find number of clusters for F and B ----
% 1)  get small set of samples --------------------------------------------
[FSampleNum,FV_num ] = size(TCData4F); % number of foreground sample
BSampleNum = size(TCData4B,1); % number of background sample

% Get Samll part of data to find number of optimal clussters --------------
FSnum_SmallPart = round(.1 *FSampleNum) ;
BSnum_SmallPart = round(.1 *BSampleNum) ;
F_SmallPartInd = round((rand(1,FSnum_SmallPart)*(FSampleNum-1))+1) ;
B_SmallPartInd = round((rand(1,BSnum_SmallPart)*(BSampleNum-1))+1) ;

TCData4F_SmallPart = TCData4F(F_SmallPartInd,:) ;
TCData4B_SmallPart = TCData4B(B_SmallPartInd,:) ;

% find number of optimum cluster for small set of data 
MaxClusterNum =4 ;
ClusterMode=1 ;
[TCnum4F,TCLBL4F,TCriteria4F] = Fun_Find_NumCluster(TCData4F_SmallPart+eps, MaxClusterNum , ClusterMode) ; % find number of clusster of F data
[TCnum4B,TCLBL4B,TCriteria4B] = Fun_Find_NumCluster(TCData4B_SmallPart+eps,MaxClusterNum , ClusterMode) ;
% -------------------------------------------------------------------------
% % 2) Use estimated  number of clussters to clussters original data --------
% EstClusterNum4F =TCnum4F ;
% EstClusterNum4B =TCnum4B ;
% ClusterMode=-1 ;
% [TCnum4F,TCLBL4F,TCriteria4F] = Fun_Find_NumCluster(TCData4F+eps, EstClusterNum4F , ClusterMode) ; % find number of clusster of F data
% [TCnum4B,TCLBL4B,TCriteria4B] = Fun_Find_NumCluster(TCData4B+eps,EstClusterNum4B , ClusterMode) ;
if FunMod==1 % Use same clusster sizes ====================================
% 3) push clussters to have same number of samples ------------------------
SSC_NumS=round(max(FSnum_SmallPart,BSnum_SmallPart)) ;             % Total number of samples
SSC_NumS4FC=round(max(FSnum_SmallPart,BSnum_SmallPart)/TCnum4F) ;  % number of samples for every classes
SSC_NumS4BC=round(max(FSnum_SmallPart,BSnum_SmallPart)/TCnum4B) ;  % number of samples for every classes

% =========================================================================
%  Make same size clussters 4 Foreground samples ------------------------- 
TCnum = TCnum4F ;
TClBL = TCLBL4F ; 
TFV = TCData4F_SmallPart ; % Feature vector of data ----------------------
TNumGSC = SSC_NumS4FC ; % number of generatin samples per class 

RGSFV= [] ; % keep feature vecotr of generated samples
RGSLBL=[] ; % Generated sapmple class lable 
TotalInC = [1 : length(TClBL) ]'  ; % total index of classes 
for i=1 : TCnum
    TIn4C = TotalInC(TClBL==i) ;
    
    TGSInd=round(rand (TNumGSC,1)* (length(TIn4C)-1)+1 ); % index of generated samples ; 
    TGSCLBL = zeros((TNumGSC),1) ; TGSCLBL(:)=i ; 
    RGSFV =[RGSFV; TFV(TGSInd,:)] ;
    RGSLBL = [RGSLBL ;TGSCLBL ] ; 
end
RGSFV4F = RGSFV ; % generated samples for background with same size clussters
RGSLBL4F =RGSLBL ;% class lable of generated samples for background with same size clussters

%  Make same size clussters 4 Background samples ------------------------- 
TCnum = TCnum4B ;
TClBL = TCLBL4B ; 
TFV = TCData4B_SmallPart ; % Feature vector of data ----------------------
TNumGSC = SSC_NumS4BC ; % number of generatin samples per class 

RGSFV= [] ; % keep feature vecotr of generated samples
RGSLBL=[] ; % Generated sapmple class lable 
TotalInC = [1 : length(TClBL) ]'  ; % total index of classes 
for i=1 : TCnum
    TIn4C = TotalInC(TClBL==i) ;
    
    TGSInd=round(rand (TNumGSC,1)* (length(TIn4C)-1)+1 ); % index of generated samples ; 
    TGSCLBL = zeros((TNumGSC),1) ; TGSCLBL(:)=i ; 
    RGSFV =[RGSFV; TFV(TGSInd,:)] ;
    RGSLBL = [RGSLBL ;TGSCLBL ] ; 
end

RGSFV4B = RGSFV ;  % generated samples for background with same size clussters
RGSLBL4B =RGSLBL ; % class labels of generated  samples for background with same size clussters
% Remove redundant information 
RGSFV=[] ; RGSLBL= [] ; TGSCLBL=[] ; TGSInd=[] ; 

TDataTrain_Ext = [RGSFV4F ; RGSFV4B];
TDataTrainLBl_Ext = [RGSLBL4F ; (RGSLBL4B+TCnum4F)] ;
%% ========================================================================
elseif FunMod==2  % use small set for LDA 
    
TDataTrain_Ext = [TCData4F_SmallPart ; TCData4B_SmallPart];
TDataTrainLBl_Ext = [TCLBL4F ; (TCLBL4B+TCnum4F)] ;    
    
end




% S2) Use LDA  ========================================================

options.Regu=1 ;
[LDA_eigvector, LDA_eigvalue, LDAelapse] = LDA(TDataTrainLBl_Ext,options,TDataTrain_Ext);

% clear redundant data ------------
TDataTrainLBl_Ext= [] ; TDataTrain_Ext=[] ;
TCData4F_Ext=[] ; TCData4B_Ext=[] ;
TCLBL4F_Ext=[]; TCLBL4B_Ext=[] ;
% -----------------------------------

LDA_InoThreshold=InfThr ;
LDA_sel=1 ; % find top fist eign vector whic keep .9 percent of informations
while sum(LDA_eigvalue(1:LDA_sel))/sum(LDA_eigvalue)  <LDA_InoThreshold
    LDA_sel=LDA_sel+1  ;
end


R_RDimLDA = FV_Vect* LDA_eigvector(:,1:LDA_sel) ; % reduced dimension by LDA

% for i=1 : LDA_sel
%     TempC=[] ;
%     TempC =  R_RDimLDA(:,i);
%     TempC=TempC - min(TempC(:)) ; TempC= TempC * (255 / max(TempC(:))) ;
%     R_RDimLDA(:,i)=TempC ;
%     TempC=[] ;
% end

R_LDA_EigenVal =LDA_eigvalue  ; % eigen values


%% Return Result ----------------------------------------------------------
RFVLDA =  reshape(R_RDimLDA , [FV_H,FV_W,LDA_sel]); % ------------------

RLDA_Eig = R_LDA_EigenVal ;  % return Eigen values -------------------------

%--------------------------------------------------------------------------

% Take three LDA Feature --------------------------------------------------
% SelDimInd=3 ;
% TCDataPCA_LDA = FV_Vect* LDA_eigvector(:,1:SelDimInd) ;
%
% TCDataPCA_LDAMat = reshape(TCDataPCA_LDA , [FV_H,FV_W,SelDimInd]);
% TCDataPCA_LDAMatScaled=[] ;
% for k=1 :SelDimInd
%     a=TCDataPCA_LDAMat(:,:,k) ;
%     a=a-min(a(:)) ; a=a*255/max(a(:)) ;
%     TCDataPCA_LDAMatScaled (:,:,k)=a ;
% end
% a=[] ;
%
% ActiveI = double(round(TCDataPCA_LDAMatScaled)) ;






% -------------------------------------------------------------------------





end

