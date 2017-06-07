function [RAlpha,RConf, RTAlpha ,RIFG,RIBG,RIWeight,RIEct ] = Fun_SharedDistMatting_V55_7_3 (IColor,ITexture, IGrad , MaskAct,  BlKNumSet,NumKg, EWC,EWT,EWCT,varargin )


%% V6 : in this version the histogram of blocks are used to find the E_c and
%% E_t of objective functions. 




% Input Variables  -----------------------------------
% EW1   : effectiveness weight for matting process
% EWS   : smoothness weight for smoothing part 
% varargin{1}  : Number of generated sampples per clusster in block 
% varargin{2}  : 0= random sampling  , 1=mean value (automatically reject Varargin{1} )



%% ========================================================================


RmaskB = (MaskAct==5) ;  % background mask
RmaskF = (MaskAct==1) ;  % foreground mask
[Ih, Iw, Cmod]= size(IColor) ;



% 1) Color Model Construction -------------------


NumLevel = 1 ;
NumBLKH=BlKNumSet(1) ; % number of block in vertical direction
NumBLKW=BlKNumSet(2) ; % number of block in horizontal direction
ColFormat= 'rgb' ;  % color format for further processing
NormalOpt=0 ; % 1 : normalize between [0 1] ; use [0 255] rage
NumGenSample4C=3 ; % number of generated sample per clusster
SamplingOption=1 ; % 0= random sampling  , 1 = mean sampling (in this case we just use mean of clusster as its generated sample)
tic
[GBlkModelB,GBlkModelF,BlkMaskAct,BlkVar ]= GlobalBlkBaseSampleGenerator_V7 (IColor,ITexture,RmaskB , RmaskF , NumLevel,NumBLKH,NumBLKW , ColFormat , NormalOpt , NumGenSample4C ,SamplingOption) ;

%[GBlkModelB2,IPmapB2,GBlkModelF2,IPmapF2,BPmap2 ]= GlobalBlkBaseSampleGenerator_V4 (IColor,RmaskB , RmaskF , NumLevel,NumBLKH,NumBLKW , ColFormat ) ;

'color sampling is completed in '
toc



TFBAct =(BlkMaskAct==1); % foreground block based color models
TFBlkModelSet =GBlkModelF(:,:,1)  ; %  keep set of generated samples for each block
TF_BlkSample =GBlkModelF(:,:,2)  ; %  keep  the blocks sample in color and texture spaces 


TBBAct = (BlkMaskAct==5) ; % background block mask
TBBlkModelSet = GBlkModelB(:,:,1); % keep set of generated samples for each block
TB_BlkSample = GBlkModelB(:,:,2);  %  keep  the blocks sample in color and texture spaces 
% % % % % % 

% 2) Matting Process ----------------------------
% NumKg  : number of kg lines 

% EW= [.5 1 1 .25 1 1 0 1 1 ] ;

% Tip : tempcolor is used instead of IColor
%load('.\Data\Data_Part1.mat')
%% NumSPClusster  indicate number of generated samples for each clusster inside the blocks 
if length(varargin)>0 
    NumSPClusster = varargin{1} ; 
else
    NumSPClusster=1 ; 
end 
if SamplingOption==1 
     NumSPClusster=1 ; 
end
tic
 
%[IFG , IBG, ITFG , ITBG,IRAlpha,IRTAlpha, IRobust,IRWeight] = GetFBSofGMLl1BasedKgLine_V20 (IColor,ITexture,IGrad,MaskAct,TBBAct,TBBlkModelSet,TFBAct,TFBlkModelSet,NumKg,NumSPClusster,EWC,EWT,EWCT,BlkVar) ;
[IFG , IBG, ITFG , ITBG,IRAlpha,IRTAlpha, IRobust,IRWeight, IREct] = GetFBSofGMLl1BasedKgLine_V55_7_3 (IColor,ITexture,IGrad,MaskAct,TBBAct,TBBlkModelSet,TB_BlkSample,TFBAct,TFBlkModelSet,TF_BlkSample,NumKg,NumSPClusster,EWC,EWT,EWCT,BlkVar) ;
toc

%load('.\Data\Data_Part1.mat')
%% 
 

GlobalAlpha= sum((IColor -  IBG).*(IFG -  IBG),3)./ sum((IFG -  IBG).^2,3) ;
GlobalAlpha (RmaskF)=1 ;GlobalAlpha (RmaskB)=0 ;
GlobalAlpha(GlobalAlpha <0)=0  ; GlobalAlpha(GlobalAlpha >1)=1  ;

% -------------------------------------------------------------------
IRAlpha(RmaskF)=1 ; IRAlpha(RmaskB)=0 ;  IRAlpha(IRAlpha<0)=0  ;  IRAlpha(IRAlpha>1)=1  ;
IRTAlpha(RmaskF)=1 ; IRTAlpha(RmaskB)=0 ;  IRTAlpha(IRAlpha<0)=0  ;  IRTAlpha(IRAlpha>1)=1  ;

IRobust (RmaskF)=1 ; IRobust (RmaskB)=1 ;
% % SRad = 15 ;
% % EWS = [1 1 1 1] ;
% R= Fun_AlphaSmoothing(IColor,IGrad,MaskAct,ITAlpha, IRAlpha , IRobust ,SRad,EWS );


%% Return Result 
RAlpha = IRAlpha ;
RTAlpha = IRTAlpha ;
RConf  = IRobust ;
RIFG = IFG ;
RIBG = IBG ;
RIWeight = IRWeight ;

RIEct = IREct ; % Return Automatic E_c and E_t ; 








% % % % % % end


%%

%% Combination of different levels result -------------------------








