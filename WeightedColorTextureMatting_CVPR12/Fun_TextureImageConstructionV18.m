function [RTImage , RTAlpha] = Fun_TextureImageConstructionV15 (I,IGray ,MaskAct ,erodLen, TextMode, NumLevel)
% In this function we construct Texte Image based on Waevelet transform and
% Texture image image is constructed based on biltaeral fitering and new
% set of features

%% in this version 8 : New approach to construct texture image is used and
%% wavelet decomposition is done on Color Image in diffferent levles and mean value of coefficient over 3x3 windows or higher
%% is computed as Feature vector , (just mean value of CA,CV,CH,CD are used in color space in different levels )


%% Compute Texture Image ----------
I= double(I) ;
IGray=double(IGray);
I=I/255 ;
IGray=IGray/255 ;
IOrg= I ;

[IOrgH , IOrgW,Cmode]= size(IOrg) ;
nrow = IOrgH ; ncol = IOrgW ;
IOrg = imresize(IOrg,[nrow ncol]); %IOrg=uint8(IOrg);
TSize =1;
I=IOrg  ; I= imresize(I,TSize);
IGrayOrg = IGray ;
IGray=imresize(IGray,TSize) ;

IColor = double(I) ;

[Ih,Iw,Cmod] = size(I) ;

% Get MaskAct -------------------
RmaskF= (MaskAct==1) ;
RmaskB= (MaskAct==5) ;
RmaskFOrg = RmaskF ;
RmaskBOrg = RmaskB ;
MaskActOrg = RmaskBOrg * 5 + RmaskFOrg ;   % get the marked F&B pixels
MaskActOrg(MaskActOrg==0)=3  ;
RmaskF=(RmaskF- imerode(RmaskFOrg,ones(erodLen))>0);
RmaskB=(RmaskB- imerode(RmaskBOrg,ones(erodLen))>0);
MaskAct = RmaskB * 5 + RmaskF ;   % get the marked F&B pixels
MaskAct(MaskAct==0)=3  ;
% -------------------------------


FWinSizeSet =[3 5 7 11  17  23  ];  % size of Fmaskes
FMSNum=1 ;  % number of different size of FMask
for FMSizeInd=1 : FMSNum%length(FWinSizeSet)
    % design Feature Mask----------------------
    FWinSize=FWinSizeSet(FMSizeInd) ; FWinRad=(FWinSize-1)/2 ;
    FVMask1 = zeros(FWinSize) ; FVMask1 (FWinRad+1,:)=1 ;FVMask1=FVMask1>0 ;
    FVMask2 = zeros(FWinSize) ; FVMask2 (:,FWinRad+1)=1 ;FVMask2=FVMask2>0 ;
    FVMask3 = eye(FWinSize)>0 ;
    FVMask4 = fliplr(eye(FWinSize))>0;
    FVMask5= ones(FWinSize)>0 ;
    FVMask6= ones(FWinSize)/(FWinSize^2) ;
    FWin = ones(FWinSize)>0 ;
    
    if TextMode==1
        CASet{1}= IGray;  % gray level features
    elseif TextMode==2
        CASet{1}= IColor; % gray level features
    else
        error ('TextMode is unknown it shoud be 1 for graylevel and  2 for color') ;
    end
    
    tic
    % NumLevel =3 ;
    CutThr=0.04 ; 
    for level =1 : NumLevel
        
        %ActI = double(imresize(CASet{level},TSize)) ; % select blured image with filter in size ((2^ (level+1)) +1) ;
        ActI = double(CASet{level}) ; % select blured image with filter in size ((2^ (level+1)) +1) ;
        %        ActI = IGray ;
        wname='haar' ;
        [CA,CH,CV,CD] = dwt2(double(ActI),wname);
        [CoefH , CoefW, Cmode] = size(CA) ;
        CA=log(abs(CA)+1) ; CH=log(abs(CH)+1) ; CV=log(abs(CV)+1) ;CD=log(abs(CD)+1) ;
        %CA=(abs(CA)) ; CH=abs(CH) ; CV=abs(CV) ;CD=abs(CD) ;
        
        % Cut th% coeef of tops -----------------------------------
        CutThr = .05 ;
        CA= FunCoefHistCut_V4 (CA,CutThr ) ;
        CH= FunCoefHistCut_V4 (CH,CutThr );
        CV= FunCoefHistCut_V4 (CV,CutThr );
        CD= FunCoefHistCut_V4 (CD,CutThr );
        % ---------------------------------------------------------
        %CA=CA/255 ; CH=CH/255 ; CV=CV/255 ; CD=CD/255 ; 
        CASet{level+1}= CA ;
        % Work On CA --------------------------------------------------------
        Mode=2 ; % mode of FV for CA,CH, CV, CD ------------------------------
        % --------------------------------------------------------------------
        ACoef=CA ;
        RM4=[] ;RM5=[] ; RMGradRGB=[] ; RM3=[] ; 
%         RMGradRGB (:,:,1) = abs(gradient(double(ACoef(:,:,1)))) ;
%         RMGradRGB (:,:,2) = abs(gradient(double(ACoef(:,:,2)))) ;
%         RMGradRGB (:,:,3) = abs(gradient(double(ACoef(:,:,3)))) ;
%         RMGradRGB =  FunCoefHistCut_V2 (RMGradRGB,CutThr ) ;
%         RMGradRGB=RMGradRGB/max(RMGradRGB(:)) ; 
%         RM3=RMGradRGB ;
%         
        RM4=stdfilt (ACoef , FVMask5) ;
        RM4 =  FunCoefHistCut_V2 (RM4,CutThr ) ;
        
        RM5=imfilter (double(ACoef) , double(FVMask6)) ;
        
        RM4= RM4/ max(RM4(:)); 
        RM5=RM5/ max(RM5(:)) ; 
        % gather all feature vector for CA
%         TFVNum3 = size(RM3,3) ;
        TFVNum4 = size(RM4,3) ;
        TFVNum5 = size(RM5,3) ;
        
        TFVNum =TFVNum4 +TFVNum5 ;
        % TFVNum =TFVNum4  ;
        TFV = zeros(CoefH,CoefW,TFVNum);
        % % %
        %TFV (:,:,1:TFVNum3)=RM3 ;
        TFV (:,:,1:TFVNum4)=RM4 ;
        TFV (:,:,TFVNum4+1 :TFVNum4+TFVNum5 )=RM5 ;
        
        FVCA = TFV  ;
        % temp -----------
        %                 FVCA(:,:,1:TFVNum1)=RM1 ; ;
        %                 FVCA(:,:,TFVNum1+1 :TFVNum1+TFVNum2 )=RM2 ;
        %                 FVCA(:,:,TFVNum1+TFVNum2+1:TFVNum1+TFVNum2+TFVNum3)=RM5 ;
        
        %-----------------
        
        %% Work On CH --------------------------------------------------------
        
        ACoef=CH ;
        
        %RM4=stdfilt (ACoef , FVMask5) ;
        RM5=imfilter (double(ACoef) , double(FVMask6)) ;
        
        RM5 =  FunCoefHistCut_V2 (RM5,CutThr ) ;
        RM5=RM5/ max(RM5(:)) ; 
        % gather all feature vector for CA
%        TFVNum4 = size(RM4,3) ;
        TFVNum5 = size(RM5,3) ;
        
        %TFVNum =TFVNum4 +TFVNum5 ;
        TFVNum =TFVNum5 ;
        TFV = zeros(CoefH,CoefW,TFVNum);
        % % %
%         TFV (:,:,1:TFVNum4)=RM4 ;
%         TFV (:,:,TFVNum4+1 :TFVNum4+TFVNum5 )=RM5 ;

         TFV (:,:,1 :TFVNum5 )=RM5 ;
        
        FVCH = TFV  ;
        % temp -----------
        %                 FVCH(:,:,1:TFVNum1)=RM1 ; ;
        %                 FVCH(:,:,TFVNum1+1 :TFVNum1+TFVNum2 )=RM2 ;
        %                 FVCH(:,:,TFVNum1+TFVNum2+1:TFVNum1+TFVNum2+TFVNum3)=RM5 ;
        
        %-----------------
        
        %% Work On CV --------------------------------------------------------
        
        ACoef=CV ;
        
            %RM4=stdfilt (ACoef , FVMask5) ;
        RM5=imfilter (double(ACoef) , double(FVMask6)) ;
        % gather all feature vector for CA
%        TFVNum4 = size(RM4,3) ;
        TFVNum5 = size(RM5,3) ;
        RM5 =  FunCoefHistCut_V2 (RM5,CutThr ) ;
        RM5=RM5/ max(RM5(:)) ; 
        
        %TFVNum =TFVNum4 +TFVNum5 ;
        TFVNum =TFVNum5 ;
        TFV = zeros(CoefH,CoefW,TFVNum);
        % % %
%         TFV (:,:,1:TFVNum4)=RM4 ;
%         TFV (:,:,TFVNum4+1 :TFVNum4+TFVNum5 )=RM5 ;

         TFV (:,:,1 :TFVNum5 )=RM5 ;
        
        FVCV = TFV  ;
        % temp -----------
        %                 FVCV(:,:,1:TFVNum1)=RM1 ; ;
        %                 FVCV(:,:,TFVNum1+1 :TFVNum1+TFVNum2 )=RM2 ;
        %                 FVCV(:,:,TFVNum1+TFVNum2+1:TFVNum1+TFVNum2+TFVNum3)=RM5 ;
        
        %-----------------
        
        %% Work On CD --------------------------------------------------------
        
        ACoef=CD ;
        
         %RM4=stdfilt (ACoef , FVMask5) ;
        RM5=imfilter (double(ACoef) , double(FVMask6)) ;
        % gather all feature vector for CA
%        TFVNum4 = size(RM4,3) ;
        TFVNum5 = size(RM5,3) ;
        
        %TFVNum =TFVNum4 +TFVNum5 ;
        TFVNum =TFVNum5 ;
        TFV = zeros(CoefH,CoefW,TFVNum);
        RM5 =  FunCoefHistCut_V2 (RM5,CutThr ) ;
        RM5=RM5/ max(RM5(:)) ; 
        
        % % %
%         TFV (:,:,1:TFVNum4)=RM4 ;
%         TFV (:,:,TFVNum4+1 :TFVNum4+TFVNum5 )=RM5 ;

         TFV (:,:,1 :TFVNum5 )=RM5 ;
        
        FVCD = (TFV)  ;
        % temp -----------
        %                 FVCD(:,:,1:TFVNum1)=RM1 ;
        %                 FVCD(:,:,TFVNum1+1 :TFVNum1+TFVNum2 )=RM2 ;
        %                 FVCD(:,:,TFVNum1+TFVNum2+1:TFVNum1+TFVNum2+TFVNum3)=RM5 ;
        
        % ---------------------------------------------------------------------
        
        FVNumA=size(FVCA,3);         FVNumH=size(FVCH,3);
        FVNumV=size(FVCV,3);         FVNumD=size(FVCD,3);
        FVNumT = FVNumA+FVNumH+FVNumV+FVNumD ;
        
        %                 FVLevel = single(zeros(CoefH,CoefW,4*FVNum));
        %                 FVLevel(:,:,1:FVNum)=round(FVCA) ;
        %                 FVLevel(:,:,FVNum+1:2*FVNum)=round(FVCH) ;
        %                 FVLevel(:,:,2*FVNum+1:3*FVNum)=round(FVCV) ;
        %                 FVLevel(:,:,3*FVNum+1:4*FVNum)=round(FVCD) ;
        
        %                 FVNum=size(FVCH,3);
        %
        FVLevel = single(zeros(CoefH,CoefW,FVNumT));
        FVLevel(:,:,1:FVNumA)=(FVCA) ;
        FVLevel(:,:,FVNumA+1:FVNumA+FVNumH)=(FVCH) ;
        FVLevel(:,:,FVNumA+FVNumH+1:FVNumA+FVNumH+FVNumV)=(FVCV) ;
        FVLevel(:,:,FVNumA+FVNumH+FVNumV+1:FVNumA+FVNumH+FVNumV+FVNumD)=(FVCD) ;
        
        
        %% save FV of leveles ---------------------------------------------------
        FVLevelSet{level,FMSizeInd}=single(FVLevel) ;
        FVLevel = [] ;
        FVCA=[] ; FVCV=[] ; FVCH=[] ; FVCD=[] ;
        
    end
end % eof FMSize=1 : 3
toc
% Celar redundant data --------------
RMCA=[] ; RMCH=[] ; RMCV=[] ; RMCD=[] ;
CASet=[] ; ActI=[] ;
%% Processing Level ----------------
% save('Test_BilateralTexture1.mat') ;
ProcessMode=2  ; % for simple closed form
% load('Test_BilateralTexture1.mat') ;

if ProcessMode==2 % PCA - > LDA-> Closed Form Matting
    %% ==================================================================
    %% Global Texture Classification Using PCA and LDA and ClosedFormmatting
    FV=[] ;TFVnum=0 ;
    PFVnum=0 ;
    for level=1 : NumLevel
        
        for FMSizeInd=1 : FMSNum
            Lind= PFVnum+1 ;
            PFVnum =PFVnum+ size( FVLevelSet{level,FMSizeInd},3) ;
            Temp1= FVLevelSet{level,FMSizeInd} ;
            a  =[] ; 
            for k=1 : size(Temp1,3)
                a= imresize(Temp1(:,:,k),[IOrgH,IOrgW]) ;
                a= a-min(a(:)) ; a= a*255/max(a(:)) ;
                Temp2(:,:,k)=a ; 
            end
            
            
            Uind = PFVnum ;
            TempFV(:,:,Lind : Uind)=Temp2 ;
            
            
        end
    end
    
    
    % -------------------------------------
    [TFV_H,TFV_W , TFVNum] = size(TempFV) ;
    MaskActLevel= imresize(MaskAct,[TFV_H,TFV_W]);
    MaskActLevel(MaskActLevel>=4.9)=5 ; MaskActLevel(MaskActLevel<=1.1)=1 ;
    MaskActLevel((MaskActLevel<5)&(MaskActLevel>1))=3 ;
    InfThr= .9 ; % keep InfThr percent of information during dimension reduction
    % PCA ---------------------
    % Normalize Data ----------
    %             TempFV_Vect = reshape(TempFV , [TFV_H*TFV_W ,TFVNum ] ) ;
    %             TempFV_Vect=zscore(double(TempFV_Vect)) ;
    %             TempFV= reshape(TempFV_Vect , [TFV_H , TFV_W ,TFVNum ]  ) ;TempFV_Vect=[] ;
    PCAdim = 0 ;
    while PCAdim<4
        TR_PCA=Fun_PCADimReduction (TempFV  , MaskActLevel , InfThr) ;
        PCAdim=size(TR_PCA,3);
        InfThr=InfThr*1.02
    end
    % Normilize Data ----------
    %             TFVNumPCA = size(TR_PCA,3) ;  TR_PCA_Vect = reshape(TR_PCA ,[TFV_H*TFV_W ,TFVNumPCA ]) ;
    %             TR_PCA_Vect = zscore(TR_PCA_Vect) ;
    %             TR_PCA = reshape(TR_PCA_Vect , [TFV_H,TFV_W ,TFVNumPCA ]) ; TR_PCA_Vect=[] ;
    LDA_FunMode=2 ;LDA_InfThr=1 ;
    [TR_LDA,TR_LDA_Eigen] =Fun_LDADimReduction (TR_PCA  , MaskActLevel , LDA_InfThr,LDA_FunMode) ;
    
    [NumDimPCA]= size(TR_LDA,3) ;
    [NumDimLDA]= size(TR_LDA,3) ;
    
    FVLDA {level,FMSizeInd}=TR_LDA ;
    FVLDA_Eigen{level,FMSizeInd} = TR_LDA_Eigen ;
    
    % clear redundent variables -------------------------------------------
    %     TR_LDA=[] ;
    %     TR_PCA=[] ; TR_LDA_Eigen=[] ;
    %     FVLevelSet=[] ;
    % ---------------------------------------------------------------------
    
    
    [FV_H ,FV_W,FV_num]=size(TR_LDA) ;
    
    Ch1=TR_LDA(:,:,1);
    Ch2=TR_LDA(:,:,2);
    Ch3=TR_LDA(:,:,3);
    
    TCDataPCA_LDAMatScaled(:,:,1)= Ch1 ;
    TCDataPCA_LDAMatScaled(:,:,2)= Ch2 ;
    TCDataPCA_LDAMatScaled(:,:,3)= Ch3 ;
    
    
    
    % S3 ) closed Form Matting ============================================
    Constant_Map = zeros(FV_H,FV_W)>0 ;
    Constant_Map(MaskActOrg==5)=1 ; Constant_Map((MaskActOrg==1))=1 ;
    Constant_Val = zeros(FV_H,FV_W)  ;
    Constant_Val((MaskActOrg==1))=1   ;
    Constant_Val((MaskActOrg==5))= 0 ;
    
    
    
    thr_alpha=[];  epsilon=[];  win_size=[]; levels_num=4;active_levels_num=1;
    IColor1=TCDataPCA_LDAMatScaled ;
    IColor1 = double(IColor1) - min(IColor1(:)) ;IColor1 = (IColor1)/ max(IColor1(:)) ;
    consts_map =zeros(IOrgH,IOrgW) ; consts_map((MaskActOrg==1))=1 ; consts_map((MaskActOrg==5))=1 ;
    consts_vals=zeros(IOrgH,IOrgW) ;consts_vals((MaskActOrg==1))=1 ; consts_vals((MaskActOrg==5))= 0 ;
    %TAlpha=solveAlphaC2F(IColor1,consts_map,consts_vals,levels_num,active_levels_num,thr_alpha,epsilon,win_size);
    TAlpha=[] ; 
    
    
    
    ActiveI = IColor1*255 ;
    
    
    TAlphaSet{level}=TAlpha ;
    TAlpha=abs(TAlpha);
    TAlpha(TAlpha>1)=1 ;
    TAlpha(TAlpha<0)=0 ;
    %     end % eof level mode
end  % enf of Processing Mode


%% Return Results ------------------------------------------------
RTImage = uint8(ActiveI)  ;
RTAlpha = TAlpha  ;

end % end of filename index




