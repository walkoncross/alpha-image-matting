function [BSampleModel,FSampleModel,RMaskAct,VarModel]= GlobalBlkBaseSampleGenerator_V7 (I, ITexture ,RmaskB , RmaskF , NumLevel,NumBLKH,NumBLKW , ColFormat,NormOpt,NumGenSample4C, SamplingOpt ) 
%% The version is designed exclusively for Experiment 10 

%% V7 is simillar to version 6 and moreover we return the color value of
%% blocks . 


%% in this function , we construct the global block based sampling models 
% in this function we construct the first level of models in each blocks
% and then for higher levels we use the N8 adjacency of block in lower
% levels and combine their models with respect to prior infromation and
% construct a model for higher levels .
%% in V3 we use improved global block base sampling and combine the
%% similirar distribution to produce more reliable models in each level . 
%   Input Variable 
%           I           : image 
%           RmaskB      : mask of background marked pixels
%           RmaskF      : mask of foreground marked pixels
%           NumLevel    : number of levels which used in expantion of
%           models in higher levels based of N8 adjucency 
%           NumBLKH     : number of blocks in hight direction 
%           NumBLKW     : number of blocks in width direction 
%           ColFormat   : color format RGB or Ycbcr
%           SamplingOpt  : 0= random sampling  , 1 = mean sampling (in this case we just use mean of clusster as its generated sample)
% 
%   OutPut Variable 
%           BSampleModel : contain the background model and prior in
%           different levels 
%           FSampleModel : contain the Foreground model and prior in
%           different levels 
%           IPixRelation : indicate a suitable model foreach pixel based on
%           lower level . 
%% ========================================================================
%% Construct the Global_Local Samples 
%%=========================================================================

% filename1=  'D:\Ehsan\University\NTU\Matlab WorkSpace\Picturs\Data base\TestData\troll.png'  ; 
% filename2= 'D:\Ehsan\University\NTU\Matlab WorkSpace\Picturs\Data base\TestData\Trimap3\troll.png' ; 
% IOrg= imread(filename1) ;  TSize = .5 ;
% I= imread(filename1) ; I= imresize(I,TSize);
% ITrimap= imread(filename2) ; ITrimap= imresize(ITrimap,TSize);
% RmaskB= (ITrimap < 3);  RmaskF= (ITrimap >250 );
[Ih,Iw,Cmod]= size(I) ; 

if strcmpi(ColFormat,'ycbcr')      
    IColor = double(rgb2ycbcr(uint8(I))) ; 
%     IColor (:,:,4)= IG ; 
    IGray = IColor(:,:,1) ; 
elseif strcmpi(ColFormat,'rgb')
    
    if NormOpt==1 
        IColor= double(I)/255 ;
    ITexture= double(ITexture)/255 ; 
%     IColor (:,:,4)= IG ; 
    IGray=rgb2gray(IColor) ; 
    ITGray = rgb2gray((ITexture)) ;    
        
    else
        IColor= double(I) ;
    ITexture= double(ITexture) ; 
%     IColor (:,:,4)= IG ; 
    IGray=double(rgb2gray(uint8(IColor))) ; 
    ITGray = double(rgb2gray(uint8(ITexture))) ;
    end
    
end
% construct variables for level 1 
% NumBLKW= 16 ; NumBLKH= 10 ; 



BlkW = Iw / NumBLKW ; 
BlkH = Ih/  NumBLKH ; 

% Model construction for level 1 ------------------------------------------
L1_BlkSmpSet4B = cell (NumBLKH,NumBLKW) ;  % contain the generated sample for each clusster of block for background
L1_BlkCAct4B = cell (NumBLKH,NumBLKW) ;  % contain the activated clusster for background

L1_BlkSmpSet4F = cell (NumBLKH,NumBLKW) ;  % contain the generated sample for each clusster of block for foreground
L1_BlkCAct4F = cell (NumBLKH,NumBLKW) ;  % contain the activated clusster for foreground

L1_BlkVarSet= cell(NumBLKH,NumBLKW) ;  % contain variance of each clusster wrt xolor and texture
%--------------------------------------------------------------------------


% contain the global model Background 

if (round(BlkW)~=BlkW) ||(round(BlkH)~=BlkH)
    error ('BLKSize is not suitable, Plz select BLK Numbers which produce integer block size')  
    
end
TBoard= zeros(Ih,Iw) ; 
ThrMask = .85*BlkW*BlkH ;
NumGenSample = NumGenSample4C ;  % number of generated sample per clusster
MaskAct = zeros(NumBLKH,NumBLKW) ;

SampleNumB = 10 ; SampleNumF=10 ; 
BlkNumInd= 0 ; 
tic
for i=1 : NumBLKH 
    for j=1 : NumBLKW

%         [ i j]
        
        BlkNumInd=BlkNumInd+1 ; 
        % index of block 
        TLCi = (i-1)*BlkH+1 ; TLCj = (j-1)*BlkW+1 ; 
        BRCi = i*BlkH       ; BRCj = j*BlkW   ;
        
        TBlkGVal =  IGray  (TLCi : BRCi , TLCj : BRCj); % block gray values 
        TBlkTGVal =  ITGray  (TLCi : BRCi , TLCj : BRCj); % block gray values 
        TBlkVal =  IColor  (TLCi : BRCi , TLCj : BRCj,:); % block gray values 
        TBlkTVal =  ITexture  (TLCi : BRCi , TLCj : BRCj,:); % block gray values 
        TMaskB  =  RmaskB (TLCi : BRCi , TLCj : BRCj); % block background pixel mask
        TMaskF  =  RmaskF (TLCi : BRCi , TLCj : BRCj); % block foreground pixel mask

        
        if ((sum(TMaskB(:))>=ThrMask)&&(sum(TMaskF(:))==0 )) % get definitly background block 
         MaskAct(i,j)=5 ; 
         TCPNumB= min(ceil(exp(std(TBlkGVal(:))*4)),2) ; 
         TTPNumB= min(ceil(exp(std(TBlkTGVal(:))*4)),2) ; 
         TPNumB = min(TCPNumB+TTPNumB,ceil(BlkH*BlkW*.03)) ; 
        % Use EM algorithm to find parameters of GMM ----------------------
        TSelB=[] ; TSelBOrg = [] ; 
        TBlkValD1= TBlkVal(:,:,1) ;TBlkValD2= TBlkVal(:,:,2) ;TBlkValD3= TBlkVal(:,:,3) ;        
        TBlkTValD1= TBlkTVal(:,:,1) ;TBlkTValD2= TBlkTVal(:,:,2) ;TBlkTValD3= TBlkTVal(:,:,3) ;        
        TSelB (:,1)= TBlkValD1(TMaskB); % get dimention 1 of background samples
        TSelB (:,2)= TBlkValD2(TMaskB); % get dimention 2 of background samples
        TSelB (:,3)= TBlkValD3(TMaskB); % get dimention 3 of background samples
        TSelB (:,4)= TBlkTValD1(TMaskB); % get dimention 1 of background samples
        TSelB (:,5)= TBlkTValD2(TMaskB); % get dimention 2 of background samples
        TSelB (:,6)= TBlkTValD3(TMaskB); % get dimention 3 of background samples
        TSelBOrg = TSelB ; 
        
        
        % Run EM algorithm 
        options = statset('MaxIter',200);
        TSelB = TSelB + (rand(size(TSelB))/100) ;
        if TPNumB>1 
        TModelB = gmdistribution.fit(TSelB,TPNumB,'Regularize', .01,'Options',options);        
        [BDLbl ,  nlogl] = cluster(TModelB,TSelB) ;
        else
        BDLbl = ones(size(TSelB,1),1) ; 
        end
        
        %-----------------------------------------------------------------
        
        
        if SamplingOpt==0
        
        TCSampleSet= zeros(1 , 6 ,TPNumB )-1 ; 
        ActCluster= zeros(1,TPNumB);
        for ki=1 : TPNumB
            TCNumSample = sum(BDLbl==ki) ; % number of sample inside the clusster
            if (TCNumSample  >= NumGenSample)
                ActCluster(1,ki)=1 ;
                BDLblMask = repmat(BDLbl==ki,[1,6]) ;
                BClusterData = reshape(TSelB((BDLblMask)) , [TCNumSample , 6] );
                  TRandID = round(rand(1,NumGenSample)*(TCNumSample-1))+1 ; 
                  TCSamples = BClusterData(TRandID,:) ;                                                          
                  TCSampleSet (:,:,ki)=TCSamples ; 
                  TCSamples=[] ; 
            end
        end
        
        elseif (SamplingOpt==1) % use mean value of clusster as generated sample. 
      
        TCSampleSet= zeros(1 , 6 ,TPNumB )-1 ; 
        TCVarSet= zeros(1 , 2 ,TPNumB )-1 ; 
        ActCluster= zeros(1,TPNumB);
        for ki=1 : TPNumB
            TCNumSample = sum(BDLbl==ki) ; % number of sample inside the clusster
            if (TCNumSample  >= NumGenSample)
                ActCluster(1,ki)=1 ;
                BDLblMask = repmat(BDLbl==ki,[1,6]) ;
                BClusterData = reshape(TSelB((BDLblMask)) , [TCNumSample , 6] );                  
                  BDataVar = var(BClusterData) ; 
                  TCSamples = mean(BClusterData) ;                                                          
                  TCSampleSet (:,:,ki)=TCSamples ; 
                  TCVarSet (1,1,ki)=sum(BDataVar(1:3)) ; % sum of color ariance
                  TCVarSet (1,2,ki)=sum(BDataVar(4:6)) ; % sum of texture variance
                  TCSamples=[] ; 
            end
        end           
   
        end  % Eof IF sample IF
        
                    
%       [TSampleB,TModelB]=BlockClustering (TBlkVal,TMaskB, TPNumB ,SampleNumB ); % get the GMM parameters 
        L1_BlkSmpSet4B {i,j,1}= TCSampleSet ;TCSampleSet=[] ; 
        L1_BlkSmpSet4B {i,j,2}= TSelBOrg ;TSelBOrg=[] ; 
        L1_BlkVarSet {i,j}= TCVarSet ;TCVarSet=[] ; 
        L1_BlkCAct4B {i,j}= ActCluster ;ActCluster=[] ; 
        
        
        end
        
        %%  evaluate foreground samples
%         if (i==2 )& (j==31)
%             a=1 ; 
%         end
        %%  ---------------------------------------------------------------
        if ((sum(TMaskF(:))>=ThrMask)  && (sum(TMaskB(:))==0 ))
            
         MaskAct(i,j) = 1 ;     
         TCPNumF= min(ceil(exp(std(TBlkGVal(:))*4)),3) ; 
         TTPNumF= min(ceil(exp(std(TBlkTGVal(:))*4)),3) ; 
         TPNumF = min(TCPNumF+TTPNumF,ceil(BlkW*BlkH*.03)) ;
        
        %% Use EM algorithm to find parameters of GMM ----------------------
        TSelF=[] ; TSelFOrg=[] ; 
        TBlkValD1= TBlkVal(:,:,1) ;TBlkValD2= TBlkVal(:,:,2) ;TBlkValD3= TBlkVal(:,:,3) ;
        TBlkTValD1= TBlkTVal(:,:,1) ;TBlkTValD2= TBlkTVal(:,:,2) ;TBlkTValD3= TBlkTVal(:,:,3) ;  
        
        TSelF (:,1)= TBlkValD1(TMaskF); % get dimention 1 of background samples
        TSelF (:,2)= TBlkValD2(TMaskF); % get dimention 2 of background samples
        TSelF (:,3)= TBlkValD3(TMaskF); % get dimention 3 of background samples
        TSelF (:,4)= TBlkTValD1(TMaskF); % get dimention 1 of background samples
        TSelF (:,5)= TBlkTValD2(TMaskF); % get dimention 2 of background samples
        TSelF (:,6)= TBlkTValD3(TMaskF); % get dimention 3 of background samples
        TSelFOrg=TSelF ; 
        
        
        options = statset('MaxIter',400);
        TSelF = TSelF + (rand(size(TSelF))/100); 
        if TPNumF>1 
        TModelF = gmdistribution.fit(TSelF,TPNumF,'Regularize', .01,'Options',options);
        [BDLbl ,  nlogl] = cluster(TModelF,TSelF) ;
        else
            BDLbl = ones(size(TSelF,1),1) ; 
        end
        %%---------------------------------------------------------------
        % Generate Sample For Each Clusster ..               
        
       if SamplingOpt==0
                   
        TCSampleSet= zeros(NumGenSample , 6 ,TPNumF )-1 ; 
        ActCluster= zeros(1,TPNumF);
        for ki=1 : TPNumF
            TCNumSample = sum(BDLbl==ki) ; % number of sample inside the clusster
            if (TCNumSample  >= NumGenSample)
                ActCluster(1,ki)=1 ;
                BDLblMask = repmat(BDLbl==ki,[1,6]) ;
                BClusterData = reshape(TSelF((BDLblMask)) , [TCNumSample , 6] );
                  TRandID = round(rand(1,NumGenSample)*(TCNumSample-1))+1 ; 
                  TCSamples = BClusterData(TRandID,:) ;                                                          
                  TCSampleSet (:,:,ki)=TCSamples ; 
                  TCSamples=[] ; 
            end
        end % eof for i 

       elseif SamplingOpt==1 
        TCSampleSet= zeros(1 , 6 ,TPNumF )-1 ; 
        TCVarSet= zeros(1 , 2 ,TPNumF )-1 ;
        ActCluster= zeros(1,TPNumF);
        for ki=1 : TPNumF
            TCNumSample = sum(BDLbl==ki) ; % number of sample inside the clusster
            if (TCNumSample  >= NumGenSample)
                ActCluster(1,ki)=1 ;
                BDLblMask = repmat(BDLbl==ki,[1,6]) ;
                BClusterData = reshape(TSelF((BDLblMask)) , [TCNumSample , 6] );
                BDataVar = var(BClusterData) ;  
                TCSamples = mean(BClusterData) ;                                                          
                TCSampleSet (:,:,ki)=TCSamples ; 
                TCVarSet (1,1,ki)=sum(BDataVar(1:3)) ; % sum of color ariance
                TCVarSet (1,2,ki)=sum(BDataVar(4:6)) ; % sum of  texture variance
                TCSamples=[] ; 
            end
        end % eof for i 

       end % Eof IF if SamplingOpt==1 
        L1_BlkSmpSet4F {i,j,1}= TCSampleSet ;TCSampleSet=[] ;
        L1_BlkSmpSet4F {i,j,2}= TSelFOrg ; TSelFOrg= [] ; 
        L1_BlkVarSet {i,j}= TCVarSet ;TCVarSet=[] ; 
        L1_BlkCAct4F {i,j}= ActCluster ;ActCluster=[] ; 
        
        end
        
                
    end
end
 

toc

%% Return Values ==========================================================

BSampleModel=L1_BlkSmpSet4B ;

FSampleModel=  L1_BlkSmpSet4F;

RMaskAct = MaskAct ; 
VarModel = L1_BlkVarSet; 




end % end of function 




function  [R]=SamplingAnalysis (I ,RmaskB,GmObjB)
% comparison between orginal samples and generated ones 



%%
[Ih, Iw,Cmod]= size(I) ;
IColor= I;  % work on RGB color space
SampleRGB_F=[] ; SampleRGB_B=[];
Find=0    ; Bind=0 ;

for i=1 : Cmod 
IColorD1=IColor(:,:,i);
SampleRGB_B(i,:)= IColorD1(RmaskB);
end

RGBHistB=[] ; 
for k=1 : Cmod
    RGBHistB =[RGBHistB, imhist(uint8(SampleRGB_B(k,:) ))'];

end

% construct new data based on samples

NB=sum(RmaskB(:)) ; 
RSB = random(GmObjB,NB);

RS_RGB_B = uint8(RSB)' ;
RS_RGBHistB = [imhist(uint8(RS_RGB_B(1,:)))'  ,imhist(uint8(RS_RGB_B(2,:)))'   ,imhist(uint8(RS_RGB_B(3,:)))'  ] ;

%%
figure (1); title('background samples histogram in RGB (Blue= orginal , Red= generated) ')
hold on
plot(RS_RGBHistB, 'r')
plot(RGBHistB, 'b')
%

R=1 ;


end %EOF SamplingAnalysis  Function









