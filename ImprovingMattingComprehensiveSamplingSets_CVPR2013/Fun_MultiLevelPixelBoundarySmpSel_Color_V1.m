function [ L1C4F,L1C4B ,L2C4F, L2C4B ]= Fun_MultiLevelPixelBoundarySmpSel_Color_V1(IColor , MaskAct,  BThickness, LSmoothness, mode, CutThr )
%% V2 : Real Pixel clusstering is Used 
%% peaks in histograms for F and B individually.
% -------------------------------------------------------------------------
%
%  Input Variables 
%
%       IColor      : Color images
%       ITexture      : Texture images
%       MaskAct :  1= foreground , 5 =background , 3= unknown regions
%       LStdDiv  : parameter to control relation between samples color and
%       texture deviation from centers with global standard deviation (Color threshold = std(Foreground colors)/ LstdDiv)
%
% ------------------------------------------------------------------------
% ------------------------------------------------------------------------
[Ih,Iw,Cmod]= size(IColor) ;
% [Tmod]= size(ITexture,3) ;
Tmod=0 ; 

IFV(:,:,1:Cmod)= IColor ; 
% IFV(:,:,Cmod+1 : Cmod+Tmod)= ITexture; 
[IFVH, IFVW,FVmod]= size(IFV) ; 

RmaskB = (MaskAct==5) ;
RmaskF = (MaskAct==1) ;
RmaskU = (MaskAct==3) ; 




RMDilaterFilter = ones(2*BThickness+1)>0 ;



Boundary_B =(RmaskB - imerode(RmaskB,RMDilaterFilter ))>0 ;
Boundary_F =(RmaskF - imerode(RmaskF,RMDilaterFilter ))>0 ;






% Clusstering The Data in 2 Levels ----------------------------------------

Ind1D= [1: 1: Ih*Iw] ; Ind1D = reshape(Ind1D,[Ih,Iw]) ;

IVec_FV = reshape(IFV , [Ih*Iw, FVmod]) ;

I_i = repmat ([1:1:Ih]', [1,Iw]);
I_j = repmat ([1:1:Iw], [Ih,1]);

% index of image ----------------------------------------------------------
IVec_i = I_i(:) ;
IVec_j = I_j(:) ;





%% Foreground Sample Selection ============================================
for FBTurn=1: 2
    L2_CArrat={} ; L1_CArrat={};
    if FBTurn==1
        SampSelMask = Boundary_F ;
    end
    
    if FBTurn==2
        SampSelMask = Boundary_B ;
    end
    
    FSet_Color=[] ;FSet_Texture=[] ; FSet_Ind=[] ;
    
    SampSelMaskExp = repmat(SampSelMask,[1,1,FVmod]);
 

    NumFSample = sum(SampSelMask(:)) ;
    
 
    
    FSet_FV= reshape(IFV(SampSelMaskExp), [ NumFSample, FVmod]) ;
    FSet_Color= FSet_FV(:,1:Cmod) ; 
    FSet_Texture= FSet_FV(:,end-Tmod+1:end) ; 
  
   
    FSet_Ind(:,1) = I_i(SampSelMask);   % i index of samples
    FSet_Ind(:,2) = I_j(SampSelMask);   % i index of samples
    FSet_Ind1D  = (1:1:sum(SampSelMask(:)))' ; % local indexes
    FSet_GlobalInd1D = Ind1D(SampSelMask) ;  % global indexes
    
    % ------------------------------------------------------------
    L0{1,1}= mean(FSet_FV) ; 
    L0{1,2}= std(FSet_FV) ; 
    L0{1,3}= FSet_FV ; 
    L0{1,4}= NumFSample ; 
    L0{1,5}= FSet_Ind ; 
    L0{1,6}= FSet_GlobalInd1D ; 
    Nfeatures=[Cmod, Tmod] ;
    
    CutNumSmp = max(CutThr* sum( SampSelMask(:)),2*(Cmod+Tmod));
    
    if strcmpi(mode,'CS') 
        
        ClusterMod='Color' ;
        R1 = GetNextLevelCluster(L0 , Nfeatures , ClusterMod ,LSmoothness(1) );        
        
        R1= InlineCutClusster(R1, CutNumSmp) ;
    
        CutNumSmp = max(CutThr* sum( SampSelMask(:)),(Cmod+Tmod));
        ClusterMod='Spatial' ;
        R2 = GetNextLevelCluster(R1 , Nfeatures , ClusterMod ,LSmoothness(2) );        
        R2= InlineCutClusster(R2, CutNumSmp) ;
     

    else
        
        error ; 
        
        
    end
    
    

    
    
 
    
    
    

    
    
    

    if FBTurn==1
        L1C4F = R1 ;
        L2C4F = R2 ;
        
    end
    
    if FBTurn==2
         L1C4B = R1 ;
         L2C4B = R2 ;
         
        
    end

end




end % End of Function



function R=InlineCutClusster(ActLC , CutNumSmp)

    NActLC = {} ; insert=0 ;    
    for i=1 : size(ActLC)
        if ActLC{i,4}>CutNumSmp
            insert=insert+1 ;
            NActLC(insert,:) =ActLC(i,:) ;
        end
    end
    R= NActLC ; 


end










function [RPeakNumB,RHist]= GetBlockPeak (I ,RmaskB)
% compute the number of peak in the histogram of I ; 



[Ih, Iw,Cmod]= size(I) ;
if Cmod>1 
    for chi=1 : Cmod
        
        [TPeakNumB,THist]= GetBlockPeak (I(:,:,chi) ,RmaskB) ; 
        RPeakNumB(chi,1)=TPeakNumB ; 
        RHist(chi,:) = THist ; 
    end    
elseif Cmod==1
    IColor = I ;


Bind = 0  ;  % number of F&B marked Pixels
PeakNumB=0 ;
%  1-1) get foreground and background marked pixels--------
for i=1 : Ih
    for j=1 : Iw
        if RmaskB(i,j)==1 % sample background
            Bind= Bind+1 ;
            SampleB(1,Bind)= IColor(i,j);
        end

    end
end

%  1-2) compute histogram of F&B pixels--------------------

if Bind >0
    SampleBHist (1,:)= imhist(uint8(SampleB(1,:))) ; % histogram of background samples
    SmoothSize=5 ;
    SampleBHistSmooth=SampleBHist ;
    %   1-2-1 ) smoothing the F&B histograms-------------------
    for iter=1 : 1
        for i=1 : 256
            lb= max(1,i-SmoothSize);
            hb= min(256,i+SmoothSize);
            SampleBHistSmooth (1,i)=mean(SampleBHistSmooth(1,lb:hb));
        end
    end
    
    % 1-3) find number of peaks--------------------------------
    
    PeakArea=5 ;
    SampleBHistPeak=SampleBHistSmooth ;
    PeakNumB=0 ;
    
    for i=1 : 256
        lb= max(1,i-PeakArea);
        hb= min(256,i+PeakArea);
        if (SampleBHistPeak (1,i)== max(SampleBHistSmooth(1,lb:hb)))&&  (SampleBHistPeak (1,i)>0)
            SampleBHistPeak(1,i)= SampleBHistSmooth(1,i) ;
            PeakNumB=PeakNumB+1 ;
        else
            SampleBHistPeak(1,i)= 0 ;
        end
    end
    
    
end % eof bind==0

%% return the result
RPeakNumB=PeakNumB ;
RHist = SampleBHist ; 
% RHist = SampleBHistSmooth ; 

end


end % eof function

function [RPeakNumB]= GetNumofHistPeak (IHist )
%% compute the number of peak in the histogram  



Bind = 0  ;  % number of F&B marked Pixels
PeakNumB=0 ;
%   
    SampleBHist = IHist ;
    SmoothSize=5 ;
    SampleBHistSmooth=SampleBHist ;
    %   1-2-1 ) smoothing the F&B histograms-------------------
    for iter=1 : 1
        for i=1 : 256
            lb= max(1,i-SmoothSize);
            hb= min(256,i+SmoothSize);
            SampleBHistSmooth (1,i)=mean(SampleBHistSmooth(1,lb:hb));
        end
    end
    
    % 1-3) find number of peaks--------------------------------    
    PeakArea=5 ;
    SampleBHistPeak=SampleBHistSmooth ;
    PeakNumB=0 ;
    
    for i=1 : 256
        lb= max(1,i-PeakArea);
        hb= min(256,i+PeakArea);
        if (SampleBHistPeak (1,i)== max(SampleBHistSmooth(1,lb:hb)))&&  (SampleBHistPeak (1,i)>0)
            SampleBHistPeak(1,i)= SampleBHistSmooth(1,i) ;
            PeakNumB=PeakNumB+1 ;
        else
            SampleBHistPeak(1,i)= 0 ;
        end
    end
    
 % eof bind==0

%% return the result
RPeakNumB=PeakNumB ;



end % e

 