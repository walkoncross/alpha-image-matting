
function [TR_FV_F , TR_FV_B,TR_FVAlpha , TRobust,TRSelWeight,TRICohenD,L2TR_FVAlpha , L2TRobust ]  =  Best_FBSelection_MLevelClusstering_V1_1 (IColor, MaskAct,  SetL24F,SetL24B, ULbl2F, ULbl2B  , EWC,SmpMod )

%%_V8 is like V7 but alpha is wighted wrt to channel distinctions between F
%%and B

%% Weights 4 Color Information ============================================
% W1=exp(-Sigma*A ) , A = ||alpha-.5|-.5| ;
% W2= exp(- (Cp - Cest) / DistFB)
% W3=\alpha * Prob2F ×+(1-\alpha)×Prob2B
% W4 =W(Spatial Block Dist)×W(Gradient Block Dist)
% W5 = exp(- |delta I| - |delta(alpha)* (F-B)| ) ;

%% Weights 4 Texture Information ==========================================
% WT1=exp(-Sigma*A ) , A = ||alpha-.5|-.5| ;
% WT2= exp(- (Cp - Cest) / DistFB)
% WT3=\alpha * Prob2F ×+(1-\alpha)×Prob2B
% WT4 = exp(- |delta I| - |delta(alpha)* (F-B)| ) ;

%% Weights 4 Color&Texture Information ====================================
% CTW1=\Calpha * Prob2TF ×+(1-\Calpha)×Prob2TB
% CTW2 : exp(- sigma * |alpha - Talpha|)
% CTW3=PlaneCTW3TVar ×PlaneCTW3CVar
%%_________________________________________________________________________
%%_________________________________________________________________________

% F-B pair with highest weight to generate F and B sample for Pixel .
%% V15 using texture Alpha To Find best samples
% Input Variable -------------------------
%            IColor         : Color Image
%            IGWeight       : Gradient Image
%            MaskAct        : indicate the activity of pixels
%                             1= Foreground , 5= Background ,
%                             3= Unknown&activated , 0 = unknown & deactive
%            BActB          : Activated background block
%            BlkModelSetB   : Distributions model for activated B blocks
%            BActF          : Activated foreground block
%            BlkModelSetF   : Distributions model for activated F blocks
%            NumKgLine      : Number of Kg lines
% ----------------------------------------
%  OutPut Variable
%            RF             : Generated Foreground Image
%            RB             : Generated Background Image
%
%% -----------------------------------------------------------
[Ih,Iw,Cmod]= size(IColor) ;
RmaskF = (MaskAct==1);
RmaskB = (MaskAct==5);
RmaskU = (MaskAct==3);

% [Tmod]= size(ITexture,3) ;
Tmod=0 ;

FVmod= Cmod+ Tmod ;
IFV = zeros(Ih,Iw,FVmod) ;
IFV(:,:,1:Cmod) = IColor ;
% IFV(:,:,end-Tmod+1:end) = ITexture ;
IFV_Vec = reshape(IFV, Ih*Iw, FVmod) ;

UFV_Vec=reshape(IFV(repmat(RmaskU,[1,1,FVmod])),[sum(RmaskU(:)), FVmod]);

I2DIndex(:,1)= reshape(repmat([1:Ih]',[1,Iw]),[Ih*Iw,1]) ;
I2DIndex(:,2)= reshape(repmat([1:Iw],[Ih,1]),[Ih*Iw,1]) ;

I2DIndU(:,1) = I2DIndex(RmaskU,1);
I2DIndU(:,2) = I2DIndex(RmaskU,2);

NumUSmp = sum(RmaskU(:)) ;
IU2DInd= zeros(Ih,Iw) ;
IU2DInd (RmaskU)= [1:1:NumUSmp] ;


% Return Planes Initialization --------------------------------------------

% REstCol4F= zeros(Ih,Iw,cmod); REstCol4B= zeros(Ih,Iw,cmod);
% REstTex4F= zeros(Ih,Iw,cmod); REstTex4B= zeros(Ih,Iw,cmod);

TRIEct= zeros(Ih,Iw,2) ;
TRSelWeight= zeros(Ih,Iw,12) ;
TR_FV_B = zeros(size(IFV)) ;
TR_FV_F = zeros(size(IFV)) ;
TR_CAlpha =  zeros(Ih,Iw) ;
TR_TAlpha =  zeros(Ih,Iw) ;
TR_FVAlpha=  zeros(Ih,Iw) ;
TRobust =  zeros(Ih,Iw) ;
TRICohenD=  zeros(Ih,Iw,2) ;


RObj = zeros(Ih,Iw) ;
RAlpha = zeros(Ih,Iw) ;
RC4F = zeros(Ih,Iw,3) ;
RC4B = zeros(Ih,Iw,3) ;
RWeights = zeros(Ih,Iw,6) ;
BInd4F= zeros(Ih,Iw)-1 ;
BInd4B= zeros(Ih,Iw)-1 ;

%% Prepare F and B clussters ----------------------------------------------
MaxLevel = length(SetL24F) ;


for l=1 : MaxLevel
    
    L24F = SetL24F{1,l} ;
    L24B = SetL24B{1,l} ;
    
    NumFSmp = size(L24F,1) ;
    NumBSmp = size(L24B,1) ;
    
    AllFSmp =[];AllFSmpInd=[] ; FCumInd=[] ; FNumCSmp=[] ;FMeanIndex=[] ;
    FCMean=[] ; FCStd=[] ;AllFSmpCInd=[] ; AllFSmpInd2D=[] ;
    FCovSet=[] ;
    FCovSet_Index=[];
    for i=1: NumFSmp
        TSmp = L24F{i,3} ;
        TSmpInd = L24F{i,6} ;
        AllFSmp = [AllFSmp;TSmp] ;
        AllFSmpInd=[AllFSmpInd ; TSmpInd ] ;
        AllFSmpInd2D=[AllFSmpInd2D ; L24F{i,5} ] ;
        
        FCovSet(:,:,i)= cov(TSmp + 2* rand(size(TSmp))) ;
        FCovSet_Index(:,:,i)= cov(L24F{i,5}+rand(size(L24F{i,5}))/2) ;
        
        FMeanIndex=[FMeanIndex; mean(L24F{i,5})];
        FNumCSmp(i) = L24F{i,4} ;
        AllFSmpCInd(end+1 : end+FNumCSmp(i))=i ; % index of samples clussters
        
        FCMean(i,:) = L24F{i,1} ;
        FCStd(i,:) = L24F{i,2} ;
        if i==1
            FCumInd(i)=0 ;
        else
            FCumInd(i)=FCumInd(i-1)+L24F{i-1,4} ;
        end
    end
    
    AllBSmp =[];AllBSmpInd=[] ; BCumInd=[] ;BNumCSmp=[] ;
    BCMean=[] ; BCStd=[] ;AllBSmpCInd=[] ; AllBSmpInd2D=[] ;BMeanIndex=[] ;
    BCovSet=[] ; BCovSet_Index=[];
    for i=1: NumBSmp
        TSmp = L24B{i,3} ;
        TSmpInd = L24B{i,6} ;
        BCMean2DInd(i,:) = mean( L24B{i,5}) ;
        AllBSmp = [AllBSmp;TSmp] ; % All B Samples
        AllBSmpInd=[AllBSmpInd ; TSmpInd ] ; % Index of All B samples
        AllBSmpInd2D=[AllBSmpInd2D ; L24B{i,5} ] ; % 2D index of samples
        BCovSet(:,:,i)= cov(TSmp + 2* rand(size(TSmp))) ;
        BCovSet_Index(:,:,i)= cov(L24B{i,5}+ rand(size(L24B{i,5}))/2);
        
        BMeanIndex=[BMeanIndex; mean(L24B{i,5})] ;
        
        BNumCSmp(i) = L24B{i,4} ;  % number of samples per clusster
        AllBSmpCInd(end+1 : end+BNumCSmp(i))=i ;
        BCMean(i,:) = L24B{i,1} ;  % Mean Value of each clusster
        BCStd(i,:) = L24B{i,2} ;   % Standard Deviation of clussters
        
        if i==1
            BCumInd(i)=0 ;
        else
            BCumInd(i)=BCumInd(i-1)+L24B{i-1,4} ;
        end
    end
    
    
    %% Use Random Sampling for unknown pixels here generate  9 different
    %% samples for each clussters
    % SetGenSmp4F {l,i}= 0
    % SetGenSmp4B {l,i} = 0
    IMask=zeros(Ih,Iw) ;
    if strcmpi(SmpMod,'mean')
        
        SetGenSmp4F {l,1} = FCMean ;
        SetGenSmp4F {l,2} = FCStd ;
        SetGenSmp4F {l,3} = 0 ;
        SetGenSmp4F {l,4} = FMeanIndex ;
        
        SetGenSmp4B {l,1} = BCMean ;
        SetGenSmp4B {l,2} = BCStd ;
        SetGenSmp4B {l,3} = 0 ;
        SetGenSmp4B {l,4} = BMeanIndex ;
        
        
    else

        for i=1 : 9
            
            % Foreground -----------------------------------------------------
            BlkFSmp=[] ; BlkFInd = [] ; BlkFCov = []; BlkFStd=[] ;BlkFCDist=[]; BlkFGradDist=[] ; BlkF1DInd=[] ;
            RandMat4F = rand(1,NumFSmp);
            RandInd = FCumInd+ ( floor(rand(1,NumFSmp).*(FNumCSmp-1))+1) ;
            RandBlkInd = AllFSmpInd(RandInd) ;  % Global index of selected pixels ---
            RandBlkInd2Di = AllFSmpInd2D(RandInd,1);
            RandBlkInd2Dj = AllFSmpInd2D(RandInd,2);            
            BlkFSmp = IFV_Vec(RandBlkInd,:) ;          BlkNFSmp=[] ;
            BlkFInd = [RandBlkInd2Di,RandBlkInd2Dj] ;
            BlkF1DInd=RandBlkInd ;
            BlkFStd=FCStd ;            
            SetGenSmp4F {l,i,1} = BlkFSmp ;
            SetGenSmp4F {l,i,2} = BlkFStd ;
            SetGenSmp4F {l,i,3} = BlkF1DInd ;
            SetGenSmp4F {l,i,4} = BlkFInd ;
            % Background ------------------------------------------
            BlkBSmp=[] ; BlkBInd = [] ; BlkBCov = []; BlkBStd=[] ;BlkBCDist=[] ; BlkBGradDist=[] ; BlkBGradDist=[] ;BlkB1DInd=[] ;
            
            RandInd = BCumInd+ ( floor(rand(1,NumBSmp).*(BNumCSmp-1))+1) ;
            RandBlkInd = AllBSmpInd(RandInd) ;
            RandBlkInd2Di = AllBSmpInd2D(RandInd,1);
            RandBlkInd2Dj = AllBSmpInd2D(RandInd,2);
            
            BlkBSmp = IFV_Vec(RandBlkInd,:) ;
            BlkBInd = [RandBlkInd2Di,RandBlkInd2Dj] ;
            BlkB1DInd=RandBlkInd ;
            BlkBStd=BCStd ;
            
            
            SetGenSmp4B {l,i,1} = BlkBSmp ;
            SetGenSmp4B {l,i,2} = BlkBStd ;
            SetGenSmp4B {l,i,3} = BlkB1DInd ;
            SetGenSmp4B {l,i,4} = BlkBInd ;

            % Show Result --------------------------------------------------
            IMask=zeros(Ih,Iw) ;
            SelMask = zeros(Ih,Iw) ;
            SelMask(BlkF1DInd)= 1 ; SelMask=SelMask>0 ;
            SelMask= imdilate(SelMask, ones(5)) ;
            
            IMask(RmaskF)=1 ; IMask(RmaskB)=1 ;
            IMask(SelMask)= 5 ;
            
            SelMask = zeros(Ih,Iw) ;
            SelMask(BlkB1DInd)= 1 ; SelMask=SelMask>0 ;
            SelMask= imdilate(SelMask, ones(5)) ;

        end
        
    end  % EOF IF SmpMod
    
end % EOF Levels ----------------------------------------------------------
% 
if strcmpi(SmpMod,'mean')
   for lf=1 : MaxLevel 
       for lb=1 : MaxLevel 
           
           TSelFSmp = SetGenSmp4F {lf,1} ; 
           TSelFStd = SetGenSmp4F {lf,2} ; 
           
           TSelBSmp = SetGenSmp4B {lb,1} ; 
           TSelBStd = SetGenSmp4B {lb,2} ; 
           
           NumGenSmp4F = size(TSelFSmp,1) ; 
           NumGenSmp4B = size(TSelBSmp,1) ; 
          
           PlaneFV_F =  reshape(repmat(TSelFSmp , [NumGenSmp4B,1]) ,[NumGenSmp4F,NumGenSmp4B,FVmod] ) ;
           PlaneFV_B =  reshape(repmat(TSelBSmp(:)' , [NumGenSmp4F,1]) ,[NumGenSmp4F,NumGenSmp4B,FVmod] ) ;
           
               Mod = 'CohenD' ;
            CohenD_Color = Fun_CohenDOverlap_V1 (TSelBSmp(:,1:Cmod),TSelBStd(:,1:Cmod) ,ones(1,NumGenSmp4B)*10,TSelFSmp(:,1:Cmod),TSelFStd(:,1:Cmod) ,ones(1,NumGenSmp4F)*10 , Mod ) ;

           SetPlaneFV_F{lf,lb} = PlaneFV_F;
           SetPlaneFV_B{lf,lb} = PlaneFV_B ; 
           SetCohenD_Color{lf,lb}=CohenD_Color ; 

       end
   end
   
    
    
end



%

a= repmat([1 2 3 ; 4 5 6 ;  7 8 9],[ceil(Ih/3) , ceil(Iw/3)])
INNInd = a(1:Ih, 1:Iw) ;


for j=1 : Iw

    for i=1 : Ih
        
        if (RmaskU(i,j)==1)
            PInd = (j-1)*Ih + i ;
          
            PNInd=INNInd(i,j);
            PLInd4F=ULbl2F(i,j);
            PLInd4B=ULbl2B(i,j);
            
            
            if strcmpi(SmpMod,'mean')
                
                PlaneFV_F = SetPlaneFV_F{PLInd4F,PLInd4B} ;
                PlaneFV_B = SetPlaneFV_B{PLInd4F,PLInd4B} ;
                CohenD_Color=SetCohenD_Color{PLInd4F,PLInd4B} ;
                
                TSelF2DInd = SetGenSmp4F{PLInd4F,4} ;
                TSelB2DInd=  SetGenSmp4B{PLInd4B,4} ;
                
                NumGenSmp4F= size(PlaneFV_F,1) ;
                NumGenSmp4B= size(PlaneFV_F,2) ;
                
            else
            
             % Fetch F and B selected samples -----------------------------
            TSelFSmp = SetGenSmp4F {PLInd4F,PNInd,1} ;
            TSelFStd = SetGenSmp4F {PLInd4F,PNInd,2};
            TSelF1DInd = SetGenSmp4F {PLInd4F,PNInd,3};
            TSelF2DInd = SetGenSmp4F {PLInd4F,PNInd,4};
            NumGenSmp4F = size(TSelFSmp,1);
            
            TSelBSmp = SetGenSmp4B {PLInd4B,PNInd,1};
            TSelBStd = SetGenSmp4B {PLInd4B,PNInd,2};
            TSelB1DInd = SetGenSmp4B {PLInd4B,PNInd,3};
            TSelB2DInd = SetGenSmp4B {PLInd4B,PNInd,4};
            NumGenSmp4B = size(TSelBSmp,1);
            
            
            PlaneFV_F =  reshape(repmat(TSelFSmp , [NumGenSmp4B,1]) ,[NumGenSmp4F,NumGenSmp4B,FVmod] ) ;
            PlaneFV_B =  reshape(repmat(TSelBSmp(:)' , [NumGenSmp4F,1]) ,[NumGenSmp4F,NumGenSmp4B,FVmod] ) ;
            
            
      
            Mod = 'CohenD' ;
           CohenD_Color = Fun_CohenDOverlap_V1 (TSelBSmp(:,1:Cmod),TSelBStd(:,1:Cmod) ,ones(1,NumGenSmp4B)*10,TSelFSmp(:,1:Cmod),TSelFStd(:,1:Cmod) ,ones(1,NumGenSmp4F)*10 , Mod ) ;
            
           
            end
            PVal = reshape(IFV(i,j,:),[1,FVmod]) ;
            PlaneFV_P =  reshape(repmat(PVal ,[NumGenSmp4F*NumGenSmp4B , 1]),[NumGenSmp4F,NumGenSmp4B,FVmod]) ;
             
           
            % 1)  Plane W1 for color -----------------------
%             PlaneC_Alpha = sum((PlaneFV_P(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).* (PlaneFV_F(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)),3)./ (sum((PlaneFV_F(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).^2,3)) ;
            PlaneFV_Alpha = sum((PlaneFV_P - PlaneFV_B).* (PlaneFV_F - PlaneFV_B),3)./ (sum((PlaneFV_F - PlaneFV_B).^2,3)) ;
%             PlaneTAlpha = sum((PlaneFV_P(:,:,end-Tmod+1:end) - PlaneFV_B(:,:,end-Tmod+1:end)).* (PlaneFV_F(:,:,end-Tmod+1:end) - PlaneFV_B(:,:,end-Tmod+1:end)),3)./ (sum((PlaneFV_F(:,:,end-Tmod+1:end) - PlaneFV_B(:,:,end-Tmod+1:end)).^2,3)) ;

            
            
            % -------------------------------------------------------------            
            %% Compute W2= exp(- (Cp - Cest) / DistFB)
            PlaneFV_AlphaExt = repmat(PlaneFV_Alpha,[1,1,FVmod]) ;
          
            PlaneFV_Pest = PlaneFV_AlphaExt.* PlaneFV_F+ (1-PlaneFV_AlphaExt).* PlaneFV_B ;
            PlaneFV_DistEst =     sqrt(sum((PlaneFV_P-PlaneFV_Pest).^2,3));
            PlaneFVW2= exp(-PlaneFV_DistEst) ;
%             PlaneFVW2= exp(-PlaneFV_DistEst./ min(PlaneFV_DistEst(:))) ;
            
                      
            %% W3=\alpha * Prob2F ×+(1-\alpha)×Prob2B  ---------
            
            % 1) W3 Color Version -----------------------------
            PlaneCDistPF = sqrt(sum((PlaneFV_P(:,:,1:Cmod) - PlaneFV_F(:,:,1:Cmod)).^2,3));
            PlaneCDistPB = sqrt(sum((PlaneFV_P(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).^2,3));
            PlaneCDistFB = sqrt(sum((PlaneFV_F(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).^2,3));
            
            TPlaneColorPR2F= PlaneCDistPB ./ (PlaneCDistPF+PlaneCDistPB );
            TPlaneColorPR2B= PlaneCDistPF ./ (PlaneCDistPF+PlaneCDistPB );
            PlaneFVW3 = TPlaneColorPR2F.*PlaneFV_Alpha +(1-PlaneFV_Alpha).*TPlaneColorPR2B ;
            
            
            % Spatial Energy ----------------------------------------------
            
            TDist2F= sqrt((TSelF2DInd(:,1)-i).^2 +(TSelF2DInd(:,2)-j).^2 ) ; 
            TDist2B= sqrt((TSelB2DInd(:,1)-i).^2 +(TSelB2DInd(:,2)-j).^2 ) ; 
            
            TDist2FPlane = reshape(repmat(TDist2F , [NumGenSmp4B,1]) ,[NumGenSmp4F,NumGenSmp4B] ) ;
            TDist2BPlane = reshape(repmat(TDist2B(:)' , [NumGenSmp4F,1]) ,[NumGenSmp4F,NumGenSmp4B] ) ;
            
            PlaneFVW4 = exp(-TDist2FPlane./ mean(TDist2FPlane(:))).*exp(-TDist2BPlane./ mean(TDist2BPlane(:))) ;
         
            
            %%==============================================================
            %% Objective Function Computation ==============================
            ValidTolerance=.2 ;            
            PlaneFVAlphaMask = ones(size(PlaneFV_Alpha)) ;
            PlaneFVAlphaMask(PlaneFV_Alpha>1+ValidTolerance)=0 ;
            PlaneFVAlphaMask(PlaneFV_Alpha<0-ValidTolerance)=0 ;
            PlaneFVAlphaMask(isnan(PlaneFV_Alpha))=0 ;
            
           
            
            %% -----------------------------------------------
%            PlaneFVW2= PlaneFVW2./ max(PlaneFVW2(PlaneFVAlphaMask>0)) ; 
%            PlaneFVW3= PlaneFVW3./ max(PlaneFVW3(PlaneFVAlphaMask>0)) ; 
%            PlaneFVW4= PlaneFVW4./ max(PlaneFVW4(PlaneFVAlphaMask>0)) ; 
%            
%          PlaneObj=((T2CRatio.*PlaneW2.*PlaneAlphaMask) + ((1-T2CRatio).*PlaneTW2.*PlaneTAlphaMask)).*  (PlaneW4.^EWC(4)) ;             
           CohenWeight = CohenD_Color./ max(CohenD_Color(:)) ;
          
           PlaneObj=(PlaneFVW2.^EWC(2)).*(PlaneFVW3.^EWC(3)).*(PlaneFVW4.^EWC(4)).*(CohenWeight.^EWC(1)) ; 
           PlaneObj=PlaneFVAlphaMask.*PlaneObj ; 
            
            % ------------------------------------------------
            if NumFSmp==1
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                SelIndOffsetF   = 1  ;               % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB   = TMaxInd  ;         % offset of selected B (btween all selected GBs) ;
            elseif NumBSmp==1
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                SelIndOffsetF = TMaxInd ;            % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB = 1          ;         % offset of selected B (btween all selected GBs) ;
            else
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                [TMax2 , TMaxInd2]= max(TMax)    ;
                SelIndOffsetF =TMaxInd(TMaxInd2)  ;  % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB = TMaxInd2          ;  % offset of selected B (btween all selected GBs) ;
            end
            
            
            
            
            if sum(sum(PlaneObj>0))>= 3
                SelFVAlpha = PlaneFV_Alpha(SelIndOffsetF,SelIndOffsetB) ;
%                 SelTAlpha = PlaneTAlpha(SelIndOffsetF,SelIndOffsetB) ;
%                 SelFVAlpha =PlaneFV_Alpha(SelIndOffsetF,SelIndOffsetB) ;
                SelRobust = (PlaneObj(SelIndOffsetF,SelIndOffsetB)) ;
                
                TRICohenD(i,j,1) = CohenD_Color(SelIndOffsetF,SelIndOffsetB) ;
                
                
                
                % SelRobust = (PlaneW9(SelIndOffsetF,SelIndOffsetB)).*PlaneW6(SelIndOffsetF,SelIndOffsetB);
            else
                SelCAlpha = 0.5 ;
                SelTAlpha=0.5 ;
                SelRobust = 0 ;
                SelFVAlpha=0.5;
            end
            
            if sum(isnan(PlaneFV_Alpha(:)))>1
                a=1 ;
            elseif sum(isnan(PlaneObj(:)))>1
                a=1 ;
            end
            if sum(SelRobust(:)<0)>1
                a=1 ;
            end
            
            
            TR_FVAlpha(i,j) = SelFVAlpha ; SelFVAlpha=[] ;
            TRobust(i,j) = SelRobust ; SelRobust=[] ;
            
            TR_FV_F(i,j,:) = PlaneFV_F(SelIndOffsetF,SelIndOffsetB,:);
            TR_FV_B(i,j,:) = PlaneFV_B(SelIndOffsetF,SelIndOffsetB,:);
            
            TRSelWeight(i,j,1) = CohenWeight(SelIndOffsetF,SelIndOffsetB) ;
            TRSelWeight(i,j,2) = PlaneFVW2(SelIndOffsetF,SelIndOffsetB) ;
            TRSelWeight(i,j,3) = PlaneFVW3(SelIndOffsetF,SelIndOffsetB) ;
            TRSelWeight(i,j,4) = PlaneFVW4(SelIndOffsetF,SelIndOffsetB) ;
            
%             TRSelWeight(i,j,3) = PlaneW3(SelIndOffsetF,SelIndOffsetB) ;
%             TRSelWeight(i,j,4) = PlaneW4(SelIndOffsetF,SelIndOffsetB) ;

            
            %% *************************************************************
            
            
        end % Enof RmaskU==1
        
    end
    j
end

TR_FVAlpha(RmaskF)=1 ; TR_FVAlpha(RmaskB)=0 ;



% Level 2 Processing ------------------------------------------------------
RmaskF_Exp = repmat(RmaskF,[1,1,3])  ; 
RmaskB_Exp = repmat(RmaskB,[1,1,3])  ; 
TR_FV_F (RmaskF_Exp) = IColor(RmaskF_Exp) ;
TR_FV_B (RmaskB_Exp) = IColor(RmaskB_Exp) ;


L2TR_FVAlpha = zeros(Ih,Iw) ; 

L2TRobust = zeros(Ih,Iw) ; 



NDist=3 ; 

for j=1 : Iw
    for i=1 : Ih
        
        
              if (RmaskU(i,j)==1)
                  
                  
                  i_1 = max(1,i-NDist ) ;              
                  i_2 = min(Ih,i+NDist) ; 
                  j_1 = max(1,j-NDist)  ;                   
                  j_2 = min(Iw,j+NDist) ; 
        
                  
                  
                  TBlkSmpF = TR_FV_F(i_1:i_2 , j_1 : j_2,:) ; 
                  TBlkSmpB = TR_FV_B(i_1:i_2 , j_1 : j_2,:) ; 
                  
                  TBlkBActMsk   = RmaskU(i_1:i_2 , j_1 : j_2) | RmaskB(i_1:i_2 , j_1 : j_2) ;
                  TBlkFActMsk   = RmaskU(i_1:i_2 , j_1 : j_2) | RmaskF(i_1:i_2 , j_1 : j_2) ;
                  
                  NumAllFSmp = size(TBlkSmpF,1)*size(TBlkSmpF,2) ;
                  NumAllBSmp = size(TBlkSmpB,1)*size(TBlkSmpB,2) ; 
                  
                  TBlkSmpB_Vec = reshape(TBlkSmpB , [NumAllBSmp, FVmod]) ; 
                  TBlkSmpF_Vec = reshape(TBlkSmpF , [NumAllFSmp, FVmod]) ; 
                  
                  NumFSmp = sum(TBlkFActMsk(:)) ; 
                  NumBSmp = sum(TBlkBActMsk(:)) ;
                  
                
                  TSelFSmp = reshape(TBlkSmpF_Vec( repmat(TBlkFActMsk(:),[1,3])),[NumFSmp,FVmod]); 
                  TSelBSmp = reshape(TBlkSmpB_Vec( repmat(TBlkBActMsk(:),[1,3])),[NumBSmp,FVmod]); 
                  
                  
                  
               PlaneFV_F =  reshape(repmat(TSelFSmp , [NumBSmp,1]) ,[NumFSmp,NumBSmp,FVmod] ) ;
               PlaneFV_B =  reshape(repmat(TSelBSmp(:)' , [NumFSmp,1]) ,[NumFSmp,NumBSmp,FVmod] ) ;
                  
                     PVal = reshape(IFV(i,j,:),[1,FVmod]) ;
            PlaneFV_P =  reshape(repmat(PVal ,[NumFSmp*NumBSmp , 1]),[NumFSmp,NumBSmp,FVmod]) ;
             
           
           
            PlaneFV_Alpha = sum((PlaneFV_P - PlaneFV_B).* (PlaneFV_F - PlaneFV_B),3)./ (sum((PlaneFV_F - PlaneFV_B).^2,3)) ;
            %% Compute W2= exp(- (Cp - Cest) / DistFB)
            PlaneFV_AlphaExt = repmat(PlaneFV_Alpha,[1,1,FVmod]) ;
          
            PlaneFV_Pest = PlaneFV_AlphaExt.* PlaneFV_F+ (1-PlaneFV_AlphaExt).* PlaneFV_B ;
            PlaneFV_DistEst =     sqrt(sum((PlaneFV_P-PlaneFV_Pest).^2,3));
%             PlaneFVW2= exp(-PlaneFV_DistEst) ;
            PlaneFVW2= exp(-PlaneFV_DistEst./ mean(PlaneFV_DistEst(:))) ;
            
                      
            %% W3=\alpha * Prob2F ×+(1-\alpha)×Prob2B  ---------
            
            % 1) W3 Color Version -----------------------------
            PlaneCDistPF = sqrt(sum((PlaneFV_P(:,:,1:Cmod) - PlaneFV_F(:,:,1:Cmod)).^2,3));
            PlaneCDistPB = sqrt(sum((PlaneFV_P(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).^2,3));
            PlaneCDistFB = sqrt(sum((PlaneFV_F(:,:,1:Cmod) - PlaneFV_B(:,:,1:Cmod)).^2,3));
            
            TPlaneColorPR2F= PlaneCDistPB ./ (PlaneCDistPF+PlaneCDistPB );
            TPlaneColorPR2B= PlaneCDistPF ./ (PlaneCDistPF+PlaneCDistPB );
            PlaneFVW3 = TPlaneColorPR2F.*PlaneFV_Alpha +(1-PlaneFV_Alpha).*TPlaneColorPR2B ;
         
            
            %%==============================================================
            %% Objective Function Computation ==============================
            ValidTolerance=.2 ;            
            PlaneFVAlphaMask = ones(size(PlaneFV_Alpha)) ;
            PlaneFVAlphaMask(PlaneFV_Alpha>1+ValidTolerance)=0 ;
            PlaneFVAlphaMask(PlaneFV_Alpha<0-ValidTolerance)=0 ;
            PlaneFVAlphaMask(isnan(PlaneFV_Alpha))=0 ;

           PlaneObj=(PlaneFVW2) ; 
           PlaneObj=PlaneFVAlphaMask.*PlaneObj ; 
           
            % ------------------------------------------------
            if NumFSmp==1
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                SelIndOffsetF   = 1  ;               % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB   = TMaxInd  ;         % offset of selected B (btween all selected GBs) ;
            elseif NumBSmp==1
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                SelIndOffsetF = TMaxInd ;            % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB = 1          ;         % offset of selected B (btween all selected GBs) ;
            else
                [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                [TMax2 , TMaxInd2]= max(TMax)    ;
                SelIndOffsetF =TMaxInd(TMaxInd2)  ;  % offset of selected F (btween all selected GFs) ;
                SelIndOffsetB = TMaxInd2          ;  % offset of selected B (btween all selected GBs) ;
            end
            
            
            
            
            if sum(sum(PlaneObj>0))>= 3
                SelFVAlpha = PlaneFV_Alpha(SelIndOffsetF,SelIndOffsetB) ;
%                 SelTAlpha = PlaneTAlpha(SelIndOffsetF,SelIndOffsetB) ;
%                 SelFVAlpha =PlaneFV_Alpha(SelIndOffsetF,SelIndOffsetB) ;
                SelRobust = (PlaneObj(SelIndOffsetF,SelIndOffsetB)) ;
  
                % SelRobust = (PlaneW9(SelIndOffsetF,SelIndOffsetB)).*PlaneW6(SelIndOffsetF,SelIndOffsetB);
            else
                SelCAlpha = 0.5 ;
                SelTAlpha=0.5 ;
                SelRobust = 0 ;
                SelFVAlpha=0.5;
            end
            
            if sum(isnan(PlaneFV_Alpha(:)))>1
                a=1 ;
            elseif sum(isnan(PlaneObj(:)))>1
                a=1 ;
            end
            if sum(SelRobust(:)<0)>1
                a=1 ;
            end
        
            
            
            
            L2TR_FVAlpha(i,j) = SelFVAlpha ; SelFVAlpha=[] ;
            L2TRobust(i,j) = SelRobust ; SelRobust=[] ;
            
            L2TR_FV_F(i,j,:) = TSelFSmp(SelIndOffsetF,:);
            L2TR_FV_B(i,j,:) = TSelBSmp(SelIndOffsetB,:);
            
            L2TRSelWeight(i,j,2) = PlaneFVW2(SelIndOffsetF,SelIndOffsetB) ;
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
              end % EOF RmasjU==3
        
        
        
    end % EOF i
j
end % EOF j



L2TR_FVAlpha(RmaskF)=1 ; L2TR_FVAlpha(RmaskB)=0 ;    
