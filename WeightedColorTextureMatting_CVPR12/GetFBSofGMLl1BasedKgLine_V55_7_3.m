function [RF , RB,RTF , RTB ,RIAlpha,RITAlpha , RIRobust,RIWight,RIEct ]  =  GetFBSofGMLl1BasedKgLine_V55_7_3 (IColor,ITexture, IGWeight,MaskAct, BActOrgB ,BlkModelSetB,B_Blk_Sample ,BActOrgF ,BlkModelSetF,F_Blk_Sample,NumKgLine,NumSamPerCluster,EWC,EWT,EWCT,BlkVar )
%% V21 :  This version is designed exclusively for Experiment 10 and it is
%% similar to V20 
%% morover here we automatically estimate E_c and E_t for objective
%% function :)
%%_________________________________________________________________________
%%_________________________________________________________________________

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
% BAct = (GBlkModelF{2,1} > 0 ) ;
% BlkModelSet = GBlkModelF {1,1} ; % keep model of blocks
BActF = BActOrgF > 0 ;
BActB = BActOrgB > 0 ;
[NumBh , NumBw]= size(BActB) ;
[Ih,Iw,Cmod] = size(IColor) ;
IColor=double(IColor) ;

% Define Image of Filter Index --------------------------------------------
% Define 9 different orientation for kgline selection
TFsample = [2 4 7 ; 5,9,1 ; 8,6,3 ] ;
TExtH = (Ih + 3- mod(Ih ,3) ) / 3 ;
TExtW = (Iw + 3- mod(Iw ,3) ) / 3 ;
TExtFilterInd= repmat(TFsample, [TExtH,TExtW]) ;
IFilterInd= TExtFilterInd (1:Ih , 1:Iw) ; % index of filter which used for selection of Distribution orientation
%------------------------------------------------


BlkW = Iw / NumBw ; % block width
BlkH = Ih/  NumBh ; % block Height

Num_GS = BlkH* BlkW ; % number of generated sample per block ;
% construct block image gradiant ------------------------------------------
% which indicate average gradient in each  block
IBGWeight = zeros(NumBh,NumBw) ;
for i=1 : NumBh
    for j=1 : NumBw
        TLCi = (i-1)*BlkH+1 ; TLCj = (j-1)*BlkW+1 ;
        BRCi = i*BlkH       ; BRCj = j*BlkW   ;

        TBlkGVal =  IGWeight (TLCi: BRCi , TLCj : BRCj) ;
        TBlkGValMean = mean (abs(TBlkGVal(:))) ;
        IBGWeight (i,j) = TBlkGValMean ;
        %--------------------------------------------------
        TBlkActVal =  MaskAct (TLCi: BRCi , TLCj : BRCj) ;
        if (sum(TBlkActVal(:)==1)==Num_GS) ||(sum(TBlkActVal(:)==5)==Num_GS)
            BlkMaskAct(i,j)=1 ;

        elseif (sum(TBlkActVal(:)==0)==Num_GS);
            BlkMaskAct(i,j)=-1 ;
        else
            BlkMaskAct(i,j)=3 ;%(sum(TBlkActVal(:)==3)/(BlkW*BlkH)) ;
        end

    end
end
% -------------------------------------------------------------------------

SelN8BlkInd = zeros( 9, NumKgLine);
SelN8BlkGrad = zeros( 9, NumKgLine);
SetSelBlockB = cell(NumBh,NumBw,2) ; %BlkGS = cell(NumBh,NumBw) ;
SetSelBlockF = cell(NumBh,NumBw,2) ; %BlkGS = cell(NumBh,NumBw) ;
% keep relation between vertical index and 2D index -----------------------

BVInd = (1:1:NumBh*NumBw);  BVInd = reshape(BVInd , [NumBh , NumBw]) ;
BPlaneInd_i = (1:1:NumBh)'; BPlaneInd_i = repmat(BPlaneInd_i,[1,NumBw] ) ;
BPlaneInd_j = (1:1:NumBw); BPlaneInd_j = repmat(BPlaneInd_j,[NumBh,1] ) ;
%--------------------------------------------------------------------------
kgstep = ceil(30/BlkH) ;
% Estimate Background Model For unknowns blocks ---------------------------
for i=1 : NumBh
    for j=1 : NumBw
        SelectedBlkInd=[] ;SelBlkGrad=[]; % cotain grad values of different kg lines ;
        SelN8BlkDist=[] ;
        SelN8BlkInd = zeros( 9, NumKgLine); SelN8BlkGrad = zeros( 9, NumKgLine);
        if BActB(i,j)==0
            % neighbor indexes --------------------------------------------
            NBi (1,1)= i               ; NBi (1,2)= j               ;
            NBi (2,1)= i               ; NBi (2,2)= min(NumBw,j+1)  ;
            NBi (3,1)= max(1,i-1 )     ; NBi (3,2)= j               ;
            NBi (4,1)= i               ; NBi (4,2)= max(1,j-1)      ;
            NBi (5,1)= min(NumBh,i+1)  ; NBi (5,2)= j               ;
            NBi (6,1)= max(1,i-1 )     ; NBi (6,2)= min(NumBw,j+1)  ;
            NBi (7,1)= max(1,i-1 )     ; NBi (7,2)= max(NumBw,j-1)  ;
            NBi (8,1)= min(NumBh,i+1)  ; NBi (8,2)= max(NumBw,j-1)  ;
            NBi (9,1)= min(NumBh,i+1)  ; NBi (9,2)= min(NumBw,j+1)  ;
            %--------------------------------------------------------------
            TempBAct = BActB ;
            for Ni=1 : 9
                Theta0 =(Ni-1)*pi /16 ;
                %                 NIndi = NBi (Ni,1) ;
                %                 NIndj = NBi (Ni,2) ;
                NIndi = i ;
                NIndj = j ;


                for kgi=1 :NumKgLine
                    ThetaKg = (kgi-1)*  ((2*pi) / NumKgLine) ;
                    Theta = Theta0+ThetaKg  ;

                    ChkHit=0 ; KgLen=0 ; KgLineInd_i =NIndi ;  KgLineInd_j =NIndj ;
                    Condi = (KgLineInd_i<=NumBh )&(KgLineInd_i>=1 );
                    Condj = (KgLineInd_i<=NumBw )&(KgLineInd_i>=1 );

                    KgLineInd_Prev_i = NIndi ; % indicate previous index of KgLineInd_i
                    KgLineInd_Prev_j = NIndj ;  % indicate previous index of KgLineInd_j
                    AcumGrad = 0 ; % accumulative gradient on Kg line
                    while ((ChkHit==0))
                        KgLen = KgLen+kgstep ;
                        KgLineInd_i =NIndi+ round(KgLen* cos(Theta)) ;
                        KgLineInd_j =NIndj+ round(KgLen* sin(Theta)) ;
                        Condi = (KgLineInd_i<=NumBh )&(KgLineInd_i>=1 );
                        Condj = (KgLineInd_j<=NumBw )&(KgLineInd_j>=1 );
                        if (Condi& Condj)==0
                            break ;
                        end
                        ChkHit = TempBAct (KgLineInd_i,KgLineInd_j) ;
                        AcumGrad =AcumGrad+ abs( IBGWeight(KgLineInd_i,KgLineInd_j) - IGWeight(KgLineInd_Prev_i,KgLineInd_Prev_j)) ;
                        %                         AcumGrad =AcumGrad+ abs(IBGWeight(KgLineInd_i,KgLineInd_j)) ;
                        KgLineInd_Prev_i =  KgLineInd_i ; % update index of KgLineInd_Prev_i
                        KgLineInd_Prev_j =  KgLineInd_j ;
                        %                        KGTrack (KgLineInd_i,KgLineInd_j)= kgi ;
                    end % eof while ;

                    if ChkHit==1
                        SelectedBlkInd =[SelectedBlkInd,BVInd(KgLineInd_i,KgLineInd_j)] ;
                        SelBlkGrad =[SelBlkGrad ,AcumGrad ];
                        SelN8BlkInd  (Ni,kgi) = BVInd(KgLineInd_i,KgLineInd_j) ;
                        SelN8BlkGrad (Ni,kgi) = AcumGrad ;
                        %                         TempBAct (KgLineInd_i,KgLineInd_j)= 0 ;
                    end

                end % eof Kg
            end % eof Ni

            if length(SelectedBlkInd)==0
                SelectedBlkInd=BVInd(BActB)' ;
                SelBlkGrad=ones(size(SelectedBlkInd)) ;
            end

            [SelectedBlkInd , m,n ]=unique(SelectedBlkInd);
            SelBlkGradUnique = SelBlkGrad(m) ;
            SelectedBlkInd2D=[] ;
            SelectedBlkInd2D (1,:)= BPlaneInd_i (SelectedBlkInd);
            SelectedBlkInd2D (2,:)= BPlaneInd_j (SelectedBlkInd);
            SelN8BlkDist = sqrt(((SelectedBlkInd2D(1,:)-i).^2) +((SelectedBlkInd2D(2,:)-j).^2 )) ;
            SetSelBlockB {i,j,1} = SelectedBlkInd ;
            SetSelBlockB {i,j,2} = SelectedBlkInd2D ;
            SetSelBlockB {i,j,3} = SelBlkGradUnique ;
            SetSelBlockB {i,j,4} = SelN8BlkInd ;
            SetSelBlockB {i,j,5} = SelN8BlkGrad;
            SetSelBlockB {i,j,6} = SelN8BlkDist; % eucledian distance of selected block to activated one


        elseif BActB(i,j)==1


            SetSelBlockB {i,j,1} = BVInd(i,j) ;
            SetSelBlockB {i,j,2} = [i;j] ;
            SetSelBlockB {i,j,3} = 1 ;

            SelN8BlkInd(:,1)=BVInd(i,j)  ;
            SelN8BlkGrad(:,1)= 1 ;
            SetSelBlockB {i,j,4} = SelN8BlkInd ;
            SetSelBlockB {i,j,5} = SelN8BlkGrad;
            SetSelBlockB {i,j,6} = 0 ;


            %generate samples for blocks -------------
            %             BlkModel = BlkModelSet{i,j} ;
            %             BlkGS{i,j} = round(random(BlkModel,Num_GS))';


        end % eof  BAct(i,j)==0



    end  % eof j
 
end % eof i


% Estimate Foreground Model For unknowns blocks ---------------------------
for i=1 : NumBh
    for j=1 : NumBw
        if (i==12) && (j==10 )
           a=1
        end
        SelectedBlkInd=[] ;SelBlkGrad=[]; SelN8BlkDist=[]; % cotain grad values of different kg lines ;
        SelN8BlkInd = zeros( 9, NumKgLine); SelN8BlkGrad = zeros( 9, NumKgLine);
        if BActF(i,j)==0
            % neighbor indexes --------------------------------------------
            NBi (1,1)= i               ; NBi (1,2)= j               ;
            NBi (2,1)= i               ; NBi (2,2)= min(NumBw,j+1)  ;
            NBi (3,1)= max(1,i-1 )     ; NBi (3,2)= j               ;
            NBi (4,1)= i               ; NBi (4,2)= max(1,j-1)      ;
            NBi (5,1)= min(NumBh,i+1)  ; NBi (5,2)= j               ;
            NBi (6,1)= max(1,i-1 )     ; NBi (6,2)= min(NumBw,j+1)  ;
            NBi (7,1)= max(1,i-1 )     ; NBi (7,2)= max(NumBw,j-1)  ;
            NBi (8,1)= min(NumBh,i+1)  ; NBi (8,2)= max(NumBw,j-1)  ;
            NBi (9,1)= min(NumBh,i+1)  ; NBi (9,2)= min(NumBw,j+1)  ;
            %--------------------------------------------------------------

            %             KGTrack= double(BAct) ;
            TempBAct = BActF ;
            for Ni=1 : 9
                Theta0 =(Ni-1)*pi /16 ;

                %                 NIndi = NBi (Ni,1) ;
                %                 NIndj = NBi (Ni,2) ;
                NIndi = i ;
                NIndj = j ;
                for kgi=1 :NumKgLine
                    ThetaKg = (kgi-1)*  ((2*pi) / NumKgLine) ;
                    Theta = Theta0+ThetaKg  ;

                    ChkHit=0 ; KgLen=0 ; KgLineInd_i =NIndi ;  KgLineInd_j =NIndj ;
                    Condi = (KgLineInd_i<=NumBh )&(KgLineInd_i>=1 );
                    Condj = (KgLineInd_i<=NumBw )&(KgLineInd_i>=1 );

                    KgLineInd_Prev_i = NIndi ; % indicate previous index of KgLineInd_i
                    KgLineInd_Prev_j = NIndj ;  % indicate previous index of KgLineInd_j
                    AcumGrad = 0 ; % accumulative gradient on Kg line
                    while ((ChkHit==0))
                        KgLen = KgLen+kgstep ;
                        KgLineInd_i =NIndi+ round(KgLen* cos(Theta)) ;
                        KgLineInd_j =NIndj+ round(KgLen* sin(Theta)) ;
                        Condi = (KgLineInd_i<=NumBh )&(KgLineInd_i>=1 );
                        Condj = (KgLineInd_j<=NumBw )&(KgLineInd_j>=1 );
                        if (Condi& Condj)==0
                            break ;
                        end
                        ChkHit = TempBAct (KgLineInd_i,KgLineInd_j) ;
                        AcumGrad =AcumGrad+ abs( IBGWeight(KgLineInd_i,KgLineInd_j) - IGWeight(KgLineInd_Prev_i,KgLineInd_Prev_j)) ;
                        %                         AcumGrad =AcumGrad+ abs(IBGWeight(KgLineInd_i,KgLineInd_j)) ;
                        KgLineInd_Prev_i =  KgLineInd_i ; % update index of KgLineInd_Prev_i
                        KgLineInd_Prev_j =  KgLineInd_j ;
                        %                        KGTrack (KgLineInd_i,KgLineInd_j)= kgi ;
                    end % eof while ;

                    if ChkHit==1
                        SelectedBlkInd =[SelectedBlkInd,BVInd(KgLineInd_i,KgLineInd_j)] ;
                        SelBlkGrad =[SelBlkGrad ,AcumGrad ];
                        SelN8BlkInd  (Ni,kgi) = BVInd(KgLineInd_i,KgLineInd_j) ;
                        SelN8BlkGrad (Ni,kgi) = AcumGrad ;
                        %                         TempBAct (KgLineInd_i,KgLineInd_j)= 0 ;
                    end

                end % eof Kg
            end % eof Ni
            if length(SelectedBlkInd)==0
                SelectedBlkInd=BVInd(BActF)' ;
                SelBlkGrad=ones(size(SelectedBlkInd)) ;
            end
            [SelectedBlkInd , m,n ]=unique(SelectedBlkInd);
            SelBlkGradUnique = SelBlkGrad(m) ;
            SelectedBlkInd2D=[] ;
            SelectedBlkInd2D (1,:)= BPlaneInd_i (SelectedBlkInd);
            SelectedBlkInd2D (2,:)= BPlaneInd_j (SelectedBlkInd);
            SelN8BlkDist = sqrt(((SelectedBlkInd2D(1,:)-i).^2) +((SelectedBlkInd2D(2,:)-j).^2 )) ;
            SetSelBlockF {i,j,1} = SelectedBlkInd ;
            SetSelBlockF {i,j,2} = SelectedBlkInd2D ;
            SetSelBlockF {i,j,3} = SelBlkGradUnique ;
            SetSelBlockF {i,j,4} = SelN8BlkInd ;
            SetSelBlockF {i,j,5} = SelN8BlkGrad;
            SetSelBlockF {i,j,6} = SelN8BlkDist; % eucledian distance of selected block to activated one


        elseif BActF(i,j)==1


            SetSelBlockF {i,j,1} = BVInd(i,j) ;
            SetSelBlockF {i,j,2} = [i;j] ;
            SetSelBlockF {i,j,3} = 1 ;

            SelN8BlkInd(:,1)=BVInd(i,j)  ;
            SelN8BlkGrad(:,1)= 1 ;
            SetSelBlockF {i,j,4} = SelN8BlkInd ;
            SetSelBlockF {i,j,5} = SelN8BlkGrad;
            SetSelBlockF {i,j,6} = 0; % eucledian distance of selected block to activated one



            %generate samples for blocks -------------
            %             BlkModel = BlkModelSet{i,j} ;
            %             BlkGS{i,j} = round(random(BlkModel,Num_GS))';


        end % eof  BAct(i,j)==0



    end  % eof j

end % eof i


%% -------------------------------------------------------------------------
%
% for i=1 : NumBh
%     for j=1 : NumBw
%         TBAct = double(BAct) ;
%
%         Selvect = SetSelBlock {i,j,1} ;
%
%         TBAct(i,j)=5 ;
%         TBAct(Selvect)=TBAct(Selvect)*2 ;
%
%         imagesc(TBAct); figure(gcf)
%         pause ;
%
%
%
%
%     end
% end

%% Generate Samples Based on Generated Blocks Indexes ----------------------

% IG= zeros(Ih,Iw,3) ;  % generated foregound image
IGenF= IColor ; IGenTF= ITexture  ;
IGenB= IColor ;IGenTB= ITexture ;

NumWeight=12 ; TRIWight =zeros(Ih,Iw,NumWeight) ;
NumEct=2 ; TRIEct=zeros(Ih,Iw,NumEct) ;
IR4FB = zeros(Ih,Iw,2) ; RTMFB(:,:,1) = (MaskAct==1) ; RTMFB(:,:,2) = (MaskAct==5) ;
IR4FB(RTMFB)=1 ;
IRAlpha= zeros(Ih,Iw) ; IRAlpha2= zeros(Ih,Iw);IRTAlpha= zeros(Ih,Iw) ; IRobust = zeros(Ih,Iw);


%
%  load('.\Data\Data_Part1.mat') ;

TW1=0 ;
for j=1 : NumBw
    for i=1 : NumBh
        
        if (i==23) && (j==15)
            a=1 ;
        end
        
        
        
        %------------------------------------------------------------------
        
        %% get block information ------------------------------------------
        if (BlkMaskAct(i,j)==3)
            
            SConf4F=[] ; SConf4B=[] ;
            % model of block --------------------------------------------------
            SelTotalDistVectB = SetSelBlockB {i,j,1} ; % all distribution ;
            SelTotalGradVectB = SetSelBlockB {i,j,3} ; % gradient of Kg lines
            SelAbsGradValB    = IBGWeight(SelTotalDistVectB);% absolute gradient value of blocks
            SelBlkIndSetB     = SetSelBlockB {i,j,4} ; % set of all distribution of N8 orientation
            SelBlkGradSetB    = SetSelBlockB {i,j,5} ; % set of all Gradient Distance of N8 orientation of distribution
            SelBlkDistB       = SetSelBlockB {i,j,6} ; % Distance of selected distributions to center
            SelBlkDistB= SelBlkDistB+.1 ;
            
            SelTotalDistVectF = SetSelBlockF {i,j,1} ; % all distribution ;
            SelTotalGradVectF = SetSelBlockF {i,j,3} ; % gradient of Kg lines
            SelAbsGradValF    = IBGWeight(SelTotalDistVectF) ;% absolute gradient value of blocks
            SelBlkIndSetF     = SetSelBlockF {i,j,4} ; % set of all distribution of N8 orientation
            SelBlkGradSetF    = SetSelBlockF {i,j,5} ; % set of all Gradient Distance of N8 orientation of distribution
            SelBlkDistF       = SetSelBlockF {i,j,6} ; % Distance of selected distributions to center
            SelBlkDistF=SelBlkDistF+.1 ;
            
            
            
            %% =======================================================================================================
            %% get interested block ----------------------------------------------
            TLCi = (i-1)*BlkH+1 ; TLCj = (j-1)*BlkW+1 ; % block index (Top left corner index) ;
            BRCi = i*BlkH       ; BRCj = j*BlkW   ;     % block index (Bottom right corner) ;
            TBlkVal =  IColor  (TLCi : BRCi , TLCj : BRCj,:); % block Color values
            TBlkFind= IFilterInd (TLCi : BRCi , TLCj : BRCj) ;
            TBlkTVal= ITexture (TLCi : BRCi , TLCj : BRCj,:) ; % Texture AlphaIGWeight
            TBlkGrad= IGWeight (TLCi : BRCi , TLCj : BRCj) ; % Gradient Block
            TBlkAct= MaskAct (TLCi : BRCi , TLCj : BRCj) ;  % activity
            %------------------------------------------------------------------
            TDistC2MU=[] ; TUPrior=[] ; TUpModel=[] ;
            %             TBlkValNewF = zeros(BlkH,BlkW,3) ;TBlkValNewB = zeros(BlkH,BlkW,3) ;
            TBlkValNewF =TBlkVal ; TBlkValNewB = TBlkVal;
            TBlkValNewTF =TBlkTVal ; TBlkValNewTB = TBlkTVal;
            
            %%-----------------------------------------------------------------
            
            SelBlkIndF = SelTotalDistVectF ;  % Activated Orientation of Distribution of F
            SelBlkIndB = SelTotalDistVectB ;  % Activated Orientation of Distribution of B
            
            SelBlkGradF = SelTotalGradVectF ; % Gradient for Activated Orientation of Distribution of F
            SelBlkGradB = SelTotalGradVectB ; % Gradient for Activated Orientation of Distribution of B
            
            SelBlkIndLenF = sum(SelBlkIndF>1);
            SelBlkIndLenB = sum(SelBlkIndB>1);
            %% Get the mean value of distributions -------------------------
            TMuF=[] ; TSigmaF=[] ; TPriorF=[] ; TGradF= [] ;TBlkDistF =[];TAbsGradF=[] ;TSigmaFMean=[] ;TBSF=0 ;
            TMuFUnchanged = [] ;TMuTAF=[] ;
            
            %             for k=1 : SelBlkIndLenF
            %                 BlkModel = BlkModelSetF{SelBlkIndF(1,k)} ;
            %                 AccComp = AccComp + BlkModel.NComponents ;
            %                 TMuF(PrevAccComp+1 : AccComp, :)=BlkModel.mu  ; % mean value of distribution
            %                 TMuTAF(PrevAccComp+1 : AccComp,1) = IBlkTAlpha(SelBlkIndF(1,k));
            %                 TBSF = BlkModel.Sigma ;
            %                 TSigmaF (:,:,PrevAccComp+1 : AccComp)=TBSF ;
            %                 TSigmaFMean (1,PrevAccComp+1 : AccComp) = (TBSF(1,1,:)+TBSF(2,2,:)+TBSF(3,3,:))/3 ; % get sigma value of covariance matrix
            %                 TPriorF (1,PrevAccComp+1:AccComp) = BlkModel.PComponents ;
            %                 TGradF (1,PrevAccComp+1:AccComp)  = SelBlkGradF (1,k) ;
            %                 TAbsGradF (1,PrevAccComp+1:AccComp)  = SelAbsGradValF (1,k) ;
            %                 TBlkDistF (1,PrevAccComp+1:AccComp)  = SelBlkDistF (1,k) ;
            %                 PrevAccComp = PrevAccComp + BlkModel.NComponents ;
            %             end
            TnumSmpF=0 ;
            %NumSamPerCluster = 1 ;
            TSF = [] ;TSFCVar=[] ; TSFTVar=[] ; 
            F_Blk_Sample_Set=[] ;
            for k=1 : SelBlkIndLenF
                BlkModel = BlkModelSetF{SelBlkIndF(1,k)} ;
                F_Blk_Sample_Set =[F_Blk_Sample_Set ; F_Blk_Sample{SelBlkIndF(1,k)} ] ; 
                BlkModelVar=BlkVar{SelBlkIndF(1,k)} ;
                for ki=1 : size(BlkModel,3)
                    TempGenSample =  BlkModel(:,:,ki);
                    TempSVar =  BlkModelVar(:,:,ki);
                    if sum(TempGenSample(:)>0)>1
                        TSF(TnumSmpF+1:TnumSmpF+NumSamPerCluster,:)= TempGenSample(1:NumSamPerCluster,:);
                        TGradF (1,TnumSmpF+1:TnumSmpF+NumSamPerCluster)  = SelBlkGradF (1,k) ;
                        TAbsGradF (1,TnumSmpF+1:TnumSmpF+NumSamPerCluster)  = SelAbsGradValF (1,k) ;
                        TBlkDistF (1,TnumSmpF+1:TnumSmpF+NumSamPerCluster)  = SelBlkDistF (1,k) ;
                        TSFCVar (TnumSmpF+1:TnumSmpF+NumSamPerCluster,1)  = TempSVar (1) ;
                        TSFTVar (TnumSmpF+1:TnumSmpF+NumSamPerCluster,1)  = TempSVar (2) ;
                        TnumSmpF=TnumSmpF+NumSamPerCluster ;
                    end
                end
            end
            
            %% Background--get the mean , sigma , prior info of activated B models
            TMuB=[] ; TSigmaB=[] ; TPriorB=[] ; TGradB= [] ;TBlkDistB=[] ; TAbsGradB=[];TSigmaBMean=[] ; TBSB=0;
            TMuBUnchanged=[] ;TMuTAB = [] ;BlkModelVar=[] ; 
            
            
            TnumSmpB=0 ;
            TSB = [] ;
            TSBCVar=[] ; TSBTVar=[] ; % color and texture variance in clussters of block for each sample
            B_Blk_Sample_Set=[] ;
            for k=1 : SelBlkIndLenB
                BlkModel = BlkModelSetB{SelBlkIndB(1,k)} ;
                B_Blk_Sample_Set= [B_Blk_Sample_Set; B_Blk_Sample{SelBlkIndB(1,k)} ]  ;
                BlkModelVar=BlkVar{SelBlkIndB(1,k)} ;
                for ki=1 : size(BlkModel,3)
                    TempGenSample = BlkModel(:,:,ki);
                    TempSVar =  BlkModelVar(:,:,ki);
                    if sum(TempGenSample(:)>0)>1
                        
                        TSB(TnumSmpB+1:TnumSmpB+NumSamPerCluster,:)= TempGenSample(1:NumSamPerCluster,:);
                        TGradB (1,TnumSmpB+1:TnumSmpB+NumSamPerCluster)  = SelBlkGradB (1,k) ;
                        TAbsGradB (1,TnumSmpB+1:TnumSmpB+NumSamPerCluster)  = SelAbsGradValB (1,k) ;
                        TBlkDistB (1,TnumSmpB+1:TnumSmpB+NumSamPerCluster)  = SelBlkDistB (1,k) ;
                        TSBCVar (TnumSmpB+1:TnumSmpB+NumSamPerCluster,1)  = TempSVar (1) ;
                        TSBTVar (TnumSmpB+1:TnumSmpB+NumSamPerCluster,1)  = TempSVar (2) ;
                        TnumSmpB=TnumSmpB+NumSamPerCluster ;
                        
                    end
                    
                end
            end
            %% ============================================================
            %% Compute Overlap of Foreground and background Distributions ---------
            NumBin = 25 ; BinMinSample_Percent=0.01 ; 
            Color_OverLap = FunOverLapDist_V0 (F_Blk_Sample_Set(:,1:3) , B_Blk_Sample_Set(:,1:3), NumBin, BinMinSample_Percent) ;
            Texture_OverLap = FunOverLapDist_V0 (F_Blk_Sample_Set(:,4:6) , B_Blk_Sample_Set(:,4:6), NumBin, BinMinSample_Percent) ;
            TSigmaOverlap = 3 ; 
            
%             Texture_OverLap=sqrt(Texture_OverLap) ; 
%             Color_OverLap=sqrt(Color_OverLap) ; 
            TOverlapPortion= Texture_OverLap / (Texture_OverLap +Color_OverLap) ; 
            COverlapPortion= Color_OverLap / (Texture_OverLap +Color_OverLap) ;
          EtSigma = 2 ; 
          EcSigma=1 ; 
          E_t = exp(- EtSigma*TOverlapPortion) ; 
           E_c = exp(- EcSigma*COverlapPortion) ; 
           
           

%             E_c = E_c / max(E_c,E_t) ; 
%             E_t = E_t / max(E_c,E_t) ; 
            
  
%             
            EWC(2)=E_c ;
            EWCT(1)=E_t ; 
            % -------------------------------------------------------------
            %% ============================================================
            %% Get samples for model F ------------------------------------
            TMuFNew=[] ;
            NumDistMuF=TnumSmpF ; NumDistMuB = TnumSmpB ;
            PlaneF= zeros(NumDistMuF,NumDistMuB,3);
            PlaneTF= zeros(NumDistMuF,NumDistMuB,3);
            for ci=1 : 3
                TS1Row = TSF(:,ci);% index from 1-3 contain color values
                TST1Row = TSF(:,3+ci);% index from 4-6 contain texture values
                TempS = repmat(TS1Row,[1,NumDistMuB])   ;
                PlaneF(:,:,ci)=TempS ; % contain selected color samples from background blocks
                TempS = repmat(TST1Row,[1,NumDistMuB])   ;
                PlaneTF(:,:,ci)=TempS ; % contain selected texture samples from background blocks
            end
            PlaneFCVar=repmat(TSFCVar,[1,NumDistMuB]) ; % color variance of block of  foreground sample
            PlaneFTVar=repmat(TSFTVar,[1,NumDistMuB]) ;% texture variance of block of  bforeground sample
            
            %% Get samples for model B ------------------------------------
            PlaneB= zeros(NumDistMuF,NumDistMuB,3);
            PlaneTB= zeros(NumDistMuF,NumDistMuB,3);
            TMuBNew=[] ;
            for ci=1 : 3
                TS1Row = TSB(:,ci);% index from 1-3 contain color values
                TST1Row = TSB(:,3+ci); % index from 4-6 contain texture values
                TempS = repmat(TS1Row,[1,NumDistMuF])   ;
                PlaneB(:,:,ci)=TempS'; % contain selected color samples from background blocks
                TempS = repmat(TST1Row,[1,NumDistMuF])   ;
                PlaneTB(:,:,ci)=TempS'; % contain selected texture samples from background blocks
            end
            PlaneBCVar=repmat(TSBCVar,[1,NumDistMuF])' ; % color variance of block of  background sample
            PlaneBTVar=repmat(TSBTVar,[1,NumDistMuF])' ; % texture variance of block of  background sample
            % --------------------------------------------------------------
            %% get  Gradient values of blocks  ---------------------------------
            PlaneGF=[] ; PlaneGB=[] ;
            PlaneGF=repmat(TGradF',[1,NumDistMuB]) ; % Plane of Gradient Distance of F
            PlaneGB=repmat(TGradB,[NumDistMuF,1]) ;  % Plane of Gradient Distance of B
            
            
            %% Computation Of Weight ======================================
            %% W6 = Get Distribution Intersection Weight ------------------
            %             PlaneW6 = Fun_ModelIntersectionWeight (TMuB, TSigmaB , TMuF , TSigmaF , PlaneF , PlaneB) ;
            
            % use cohen d method for overlapping of distributions ----
            %PlaneW6 = Fun_ModelIntersectionWeight_V2 (TMuB, TSigmaB , TMuF , TSigmaF , PlaneF , PlaneB) ;
            %             PlaneW6=PlaneW6./ max(PlaneW6(:)) ;
            
            %% PlaneWDist = Weight of Distance (spatial distance ) ----------------------------
            PlaneBlkDistB = repmat(TBlkDistB,[NumDistMuF,1]); % spatial distance of btween selected block and interestd block
            PlaneBlkDistWB= exp(-PlaneBlkDistB./ mean(PlaneBlkDistB(:))) ;
            PlaneBlkDistWB = PlaneBlkDistWB./ max(PlaneBlkDistWB(:));
            
            PlaneBlkDistF = repmat(TBlkDistF',[1,NumDistMuB]);
            PlaneBlkDistWF= exp(-PlaneBlkDistF./ mean(PlaneBlkDistF(:))) ;
            PlaneBlkDistWF=PlaneBlkDistWF./max(PlaneBlkDistWF(:));
            
            
            PlaneWDist = (PlaneBlkDistWB.*PlaneBlkDistWF );
            
            % PlaneWGDist = gradient Weight of Distance (gradient distance )
            PlaneBlkDistWGF = exp(-PlaneGF./ mean(PlaneGF(:))) ;PlaneBlkDistWGF =PlaneBlkDistWGF ./max(PlaneBlkDistWGF(:) ) ;
            PlaneBlkDistWGB = exp(-PlaneGB./ mean(PlaneGB(:))) ;PlaneBlkDistWGB=PlaneBlkDistWGB./max(PlaneBlkDistWGB(:)) ;
            PlaneWGDist=PlaneBlkDistWGF.*PlaneBlkDistWGB ;
            %--------------------------------------------------------------
            %
            
            %--------------------------------------------------------------
            TBlkConf4F=zeros(BlkH,BlkW) ; TBlkConf4B=zeros(BlkH,BlkW) ;
            TBlkRobust=zeros(BlkH,BlkW) ; TBlkAlpha=zeros(BlkH,BlkW) ;TBlkTAlpha=zeros(BlkH,BlkW) ;
            TBlkSelWeight=zeros(BlkH,BlkW,NumWeight) ;
            % padding of color values for selected block (for W8)---------
            Ext_i =1 ; Ext_j =1 ;
            TBlkVal_Ext = padarray(TBlkVal,[Ext_i Ext_j],'replicate','both') ;
            TBlkTVal_Ext = padarray(TBlkTVal,[Ext_i Ext_j],'replicate','both') ;
            % ---------
            for kj=1 : BlkW
                for ki =1 : BlkH
                    if (ki==4) && (kj==17)
                        a=1 ;
                    end
                    
                    
                    if TBlkAct(ki,kj)==3 % select unknown but active pixels
                        %                         [ki,kj]
                        % Construct PlaneC ---------------
                        PlaneC=[];PlaneT=[];
                        PCVAL(1) = TBlkVal(ki,kj,1) ;PCVAL(2) = TBlkVal(ki,kj,2) ;PCVAL(3) = TBlkVal(ki,kj,3) ;
                        PlaneC(:,:,1) = repmat(PCVAL(1),[NumDistMuF,NumDistMuB]) ; % color value of interested pixel
                        PlaneC(:,:,2) = repmat(PCVAL(2),[NumDistMuF,NumDistMuB]) ;
                        PlaneC(:,:,3) = repmat(PCVAL(3),[NumDistMuF,NumDistMuB]) ;
                        PTVAL(1) = TBlkTVal(ki,kj,1) ;PTVAL(2) = TBlkTVal(ki,kj,2) ;PTVAL(3) = TBlkTVal(ki,kj,3) ;
                        PlaneT(:,:,1) = repmat(PTVAL(1),[NumDistMuF,NumDistMuB]) ; % texture value of interested pixel
                        PlaneT(:,:,2) = repmat(PTVAL(2),[NumDistMuF,NumDistMuB]) ;
                        PlaneT(:,:,3) = repmat(PTVAL(3),[NumDistMuF,NumDistMuB]) ;
                        
                        
                        
                        AlphaZeroFlag=1 ; AZeroIter=1 ;
                        %% ************************************************
                        
                        %% generate sample for each distribution ----------------------
                        
                        %                         TMUDiffCF = (repmat(PCVAL,[NumDistMuF 1])-TMuF  ) ;
                        %                         TMUDistCF = sqrt(sum((TMuF - repmat(PCVAL,[NumDistMuF 1])).^2,2 ));
                        %
                        %                         TMUDiffCB = (repmat(PCVAL,[NumDistMuB 1])-TMuB );
                        %                         TMUDistCB = sqrt(sum((TMuB - repmat(PCVAL,[NumDistMuB 1])).^2,2));
                        
                        
                        %% W1)  W1=exp(-Sigma*A ) , A = ||alpha-.5|-.5| ;
                        
                        % 1)  Plane W1 for color -----------------------
                        PlaneAlpha = sum((PlaneC - PlaneB).* (PlaneF - PlaneB),3)./ (sum((PlaneF - PlaneB).^2,3)) ;
                        %                         PlaneAlpha(isnan(PlaneAlpha))=-1 ;
                        SigmaW1 = 1 ;
                        PlaneW1 = exp(-SigmaW1*abs(abs(PlaneAlpha - .5)-.5));
                        Tolerance= 0.1 ;
                        PlaneW1(PlaneAlpha>1+Tolerance)=0 ; PlaneW1(PlaneAlpha< -Tolerance)=0 ;
                        PlaneW1(isnan(PlaneAlpha))=0 ;
                        % 1)  Plane TW1 for Texture -----------------------
                        PlaneTAlpha = sum((PlaneT - PlaneTB).* (PlaneTF - PlaneTB),3)./ (sum((PlaneTF - PlaneTB).^2,3)) ;
                        PlaneTW1 = exp(-SigmaW1*abs(abs(PlaneTAlpha - .5)-.5));
                        PlaneTW1(PlaneTAlpha>1+Tolerance)=0 ; PlaneTW1(PlaneTAlpha< -Tolerance)=0 ;
                        
                        % -------------------------------------------------
                        %% Compute W2= exp(- (Cp - Cest) / DistFB)
                        %1)  W2 for color version -----------------------
                        PlaneDistCF = sqrt(sum((PlaneC-PlaneF).^2,3));
                        PlaneDistCB = sqrt(sum((PlaneC-PlaneB).^2,3));
                        PlaneDistFB = sqrt(sum((PlaneF-PlaneB).^2,3));
                        
                        PlaneAlpha3ch = repmat(PlaneAlpha,[1,1,3]) ;
                        PlaneCest = PlaneAlpha3ch.* PlaneF + (1-PlaneAlpha3ch).* PlaneB ;
                        PlaceDistEst =     sqrt(sum((PlaneC-PlaneCest).^2,3));                
                        MeanDfb  = mean(PlaneDistFB(:)) ;
                        %PlaneW2 = exp(- (sqrt(sum((PlaneC-PlaneCest).^2,3)))./(MeanDfb)) ;
                        PlaneW2 = exp(- (PlaceDistEst)./mean(PlaceDistEst(:))) ;

                        %2)  W2 for Texture version ---------------------

                        PlaneDistTFB = sqrt(sum((PlaneTF-PlaneTB).^2,3));
                        
                        PlaneTAlpha3ch = repmat(PlaneTAlpha,[1,1,3]) ;
                        PlaneTest = PlaneTAlpha3ch.* PlaneTF + (1-PlaneTAlpha3ch).* PlaneTB ;
                        PlaceDistEstT =    sqrt(sum((PlaneT-PlaneTest).^2,3));                                     
                        MeanDTfb  = mean(PlaneDistTFB(:)) ;
                        PlaneTW2 = exp(- PlaceDistEstT./mean(PlaceDistEstT(:))) ;
                       % -------------------------------------------------- 
                       %% W3=\alpha * Prob2F ×+(1-\alpha)×Prob2B  ---------
                       
                       % 1) W3 Color Version -----------------------------
                        TPlaneColorPR2F= PlaneDistCB ./ (PlaneDistCF+PlaneDistCB );
                        TPlaneColorPR2B= PlaneDistCF ./ (PlaneDistCF+PlaneDistCB );
                        PlaneW3 = TPlaneColorPR2F.*PlaneAlpha +(1-PlaneAlpha).*TPlaneColorPR2B ;
                        % 2 ) W3 Texture Version --------------------------
                        PlaneDistTF = sqrt(sum((PlaneT-PlaneTF).^2,3));
                        PlaneDistTB = sqrt(sum((PlaneT-PlaneTB).^2,3));
                        TPlaneTexturePR2F= PlaneDistTB ./ (PlaneDistTF+PlaneDistTB );
                        TPlaneTexturePR2B= PlaneDistTF ./ (PlaneDistTF+PlaneDistTB );
                        PlaneTW3 = TPlaneTexturePR2F.*PlaneTAlpha +(1-PlaneTAlpha).*TPlaneTexturePR2B ;
                        % -------------------------------------------------
                        %% W4 =W(Spatial Block Dist)×W(Gradient Block Dist)
                        % indicate the transition energy required for foreground and background samples to reach interested pixel 
                        PlaneW4 = PlaneWDist.*PlaneWGDist ; 
                        
                        %%  W5 = exp(- |delta I| - |delta(alpha)* (F-B)| ) ;
                        
                        % 1) W5 Color Version -----------------------------
                        Num_Neighbor=4 ;
                        
                        NPCVAL(1,1:3)=TBlkVal_Ext(ki+Ext_i,kj+Ext_j+1,:) ;% read color values of neighborhoods --
                        NPCVAL(2,1:3)=TBlkVal_Ext(ki+Ext_i-1,kj+Ext_j,:) ;
                        NPCVAL(3,1:3)=TBlkVal_Ext(ki+Ext_i,kj+Ext_j-1,:) ;
                        NPCVAL(4,1:3)=TBlkVal_Ext(ki+Ext_i+1,kj+Ext_j,:) ;                       
                        
                        TDeltaBlkI = sqrt(sum((NPCVAL - repmat(PCVAL,[Num_Neighbor,1])).^2,2));
                        TDiffBlkI = (NPCVAL - repmat(PCVAL,[Num_Neighbor,1])) ;
                        TempW = zeros( NumDistMuF,NumDistMuB,Num_Neighbor) ;
                        PlaneC_Ni=zeros( NumDistMuF,NumDistMuB,3 );
                        for ni=1 : Num_Neighbor
                            
                            PlaneC_Ni(:,:,1)=PlaneC(:,:,1) + TDiffBlkI(ni,1);
                            PlaneC_Ni(:,:,2)=PlaneC(:,:,2) + TDiffBlkI(ni,2);
                            PlaneC_Ni(:,:,3)=PlaneC(:,:,3) + TDiffBlkI(ni,3);                           
                            PlaneAlphaNi = sum((PlaneC_Ni - PlaneB).* (PlaneF - PlaneB),3)./ (sum((PlaneF - PlaneB).^2,3)) ;
                            TDeltaAlpha = PlaneAlphaNi - PlaneAlpha ;                            
                            TRightLeftDiff =abs(abs(TDeltaAlpha .* PlaneDistFB)-abs(TDeltaBlkI(ni))) ;
                            TRightLeftDiff=TRightLeftDiff+.01 ;
                            TempW(:,:,ni) = exp(-TRightLeftDiff./ mean(TRightLeftDiff(~isnan(PlaneAlphaNi)))) ;
                            
                        end
                        PlaneW5 = sum(TempW,3)/Num_Neighbor ;
%                         
                        % 1) W5 Texture Version -----------------------------                        
                        NPTVAL(1,1:3)=TBlkTVal_Ext(ki+Ext_i,kj+Ext_j+1,:) ;% read Texture values of neighborhoods --
                        NPTVAL(2,1:3)=TBlkTVal_Ext(ki+Ext_i-1,kj+Ext_j,:) ;
                        NPTVAL(3,1:3)=TBlkTVal_Ext(ki+Ext_i,kj+Ext_j-1,:) ;
                        NPTVAL(4,1:3)=TBlkTVal_Ext(ki+Ext_i+1,kj+Ext_j,:) ;
                        
                        TDeltaBlkTI = sqrt(sum((NPTVAL - repmat(PTVAL,[Num_Neighbor,1])).^2,2));
                        TDiffBlkTI = (NPTVAL - repmat(PTVAL,[Num_Neighbor,1])) ;
                        TempTW = zeros( NumDistMuF,NumDistMuB,Num_Neighbor) ;
                        PlaneT_Ni=zeros( NumDistMuF,NumDistMuB,3 );
                        for ni=1 : Num_Neighbor
                            
                            PlaneT_Ni(:,:,1)=PlaneT(:,:,1) + TDiffBlkTI(ni,1);
                            PlaneT_Ni(:,:,2)=PlaneT(:,:,2) + TDiffBlkTI(ni,2);
                            PlaneT_Ni(:,:,3)=PlaneT(:,:,3) + TDiffBlkTI(ni,3);
                            
                            PlaneTAlphaNi = sum((PlaneT_Ni - PlaneTB).* (PlaneTF - PlaneTB),3)./ (sum((PlaneTF - PlaneTB).^2,3)) ;                            
                            TDeltaTAlpha = PlaneTAlphaNi - PlaneTAlpha ;                            
                            TRightLeftDiff =abs(abs(TDeltaTAlpha .* PlaneDistTFB)-abs(TDeltaBlkTI(ni))) ;
                            TRightLeftDiff=TRightLeftDiff+.01 ;
                            TempTW(:,:,ni) = exp(-TRightLeftDiff./ mean(TRightLeftDiff(~isnan(PlaneTAlphaNi)))) ;
                            
                        end
                        PlaneTW5 = sum(TempTW,3)/Num_Neighbor ;
                        
                        % -------------------------------------------------

                        
                      
                        %% ================================================
                        %% Objective Function Computation------------------
                        
                        % This Mask Remove invalid alpha from further considaration 
                        PlaneAlphaMask = ones(size(PlaneAlpha)) ; 
                        PlaneAlphaMask(PlaneAlpha>1+Tolerance)=0 ;
                        PlaneAlphaMask(PlaneAlpha<0-Tolerance)=0 ;

                        % PreProcessing Of Alpha -------------------------
                        if sum(sum(PlaneW1>0))<=5
                            PlaneW1(:) = 1 ;
                            PlaneAlpha= TPlaneColorPR2F ;
                            
                        end
                        
                        %% Texture and color weights ======================
                        
                       %% CTW1=\Calpha * Prob2TF ×+(1-\Calpha)×Prob2TB  ---------                                             
                        TPlaneTexturePR2F= PlaneDistTB ./ (PlaneDistTF+PlaneDistTB );
                        TPlaneTexturePR2B= PlaneDistTF ./ (PlaneDistTF+PlaneDistTB );
                        PlaneCTW1 = TPlaneTexturePR2F.*PlaneAlpha +(1-PlaneAlpha).*TPlaneTexturePR2B ;
                        
                       %% CTW2 : exp(- sigma * |alpha - Talpha|)
                        WCT2_Sigma = 3 ; 
                        PlaneCTW2= exp(-WCT2_Sigma*abs(PlaneAlpha - PlaneTAlpha)) ;
                        
                        %% CTW3=PlaneCTW3TVar ×PlaneCTW3CVar
                       PlaneCTW3CVar= exp(-PlaneFCVar/ mean(PlaneFCVar(:))) .* exp(-PlaneBCVar/ mean(PlaneBCVar(:))) ;
                       PlaneCTW3TVar= exp(-PlaneFTVar/ mean(PlaneFTVar(:))) .* exp(-PlaneBTVar/ mean(PlaneBTVar(:))) ;
                       
                       PlaneCTW3TVar=sqrt(PlaneCTW3TVar); 
                       PlaneCTW3CVar=sqrt(PlaneCTW3CVar) ; 
                       PlaneCTW3TVar=PlaneCTW3TVar./max(PlaneCTW3TVar(:)) ; 
                       PlaneCTW3CVar=PlaneCTW3CVar./max(PlaneCTW3CVar(:)) ; 
                       
                       PlaneCTW3 = PlaneCTW3TVar.*PlaneCTW3CVar ;
                     
                                         
                        
                        
                                            
                        
                        
                        
                        
                        %% Compute Objective Function_______________________
                        % 1) objective of color --------------------------
                        e1= EWC(1) ; e2= EWC(2) ; e3= EWC(3) ; e4= EWC(4) ; e5= EWC(5) ;
                                                
                        ObjPartC1 = (PlaneW1.^e1).*(PlaneW2.^e2);
                        ObjPartC2 = (PlaneW3.^e3).*(PlaneW4.^e4) ;
                        ObjPartC3 = (PlaneW5.^e5) ;
                        PlaneObjC = ObjPartC1.*ObjPartC2.*ObjPartC3 ; 
                        PlaneObjC=PlaneObjC./max(PlaneObjC(:)) ;
                        
                        % OBjective of texture ---------------------------
                        e1= EWT(1) ; e2= EWT(2) ;  e3= EWT(3) ; e4= EWT(4) ;  
                                                
                        ObjPartT1 = (PlaneTW1.^e1).*(PlaneTW2.^e2);
                        ObjPartT2 = (PlaneTW3.^e3).*(PlaneTW5.^e4) ;                        
                        PlaneObjT = ObjPartT1.*ObjPartT2 ; 
                        PlaneObjT=PlaneObjT./max(PlaneObjT(:)) ;
                        % OBjective of texture ---------------------------                        
                        e1= EWCT(1) ; e2= EWCT(2) ;   e3= EWCT(3) ;                                              
                        ObjPartCT1 = (PlaneCTW1.^e1).*(PlaneCTW2.^e2).*(PlaneCTW3.^e3);
                        PlaneObjCT = ObjPartCT1 ; 
                        PlaneObjCT=PlaneObjCT./ max(PlaneObjCT(:)) ; 
                        %% ------------------------------------------------
                        PlaneObj= PlaneObjC.*PlaneObjT.*PlaneObjCT ; 
                        PlaneObj = PlaneObj.*PlaneAlphaMask ; 
                        
                        
                        %% finde index of similar F&B distribution --------
                        if NumDistMuF==1
                            [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                            SelIndOffsetF   = 1  ;               % offset of selected F (btween all selected GFs) ;
                            SelIndOffsetB   = TMaxInd  ;         % offset of selected B (btween all selected GBs) ;
                        elseif NumDistMuB==1
                            [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                            SelIndOffsetF = TMaxInd ;            % offset of selected F (btween all selected GFs) ;
                            SelIndOffsetB = 1          ;         % offset of selected B (btween all selected GBs) ;
                        else
                            [TMax , TMaxInd] = max(PlaneObj) ;   % find  indexes i,j of maximum value
                            [TMax2 , TMaxInd2]= max(TMax)    ;
                            SelIndOffsetF =TMaxInd(TMaxInd2)  ;  % offset of selected F (btween all selected GFs) ;
                            SelIndOffsetB = TMaxInd2          ;  % offset of selected B (btween all selected GBs) ;
                            
                        end
                        %--------------------------------------------------------------
                        % get the mean and sigma of selected distribution
                        SelF  = TSF(SelIndOffsetF,1:3) ;
                        SelTF  = TSF(SelIndOffsetF,4:6) ;
                        
                        SelB  = TSB(SelIndOffsetB,1:3) ;
                        SelTB  = TSB(SelIndOffsetB,4:6) ;
                                               
                        
                        
                        if sum(sum(PlaneObj>0))>= 3
                            SelAlpha = PlaneAlpha(SelIndOffsetF,SelIndOffsetB) ;
                            SelTAlpha = PlaneTAlpha(SelIndOffsetF,SelIndOffsetB) ;
                            %SelRobust =  (PlaneW6(SelIndOffsetF,SelIndOffsetB)*PlaneW7(SelIndOffsetF,SelIndOffsetB)).^2;
                            %SelRobust =PlaneW6(SelIndOffsetF,SelIndOffsetB).*PlaneW8(SelIndOffsetF,SelIndOffsetB) ;
                            SI4F = SelIndOffsetF ; SI4B = SelIndOffsetB ;
                            SelRobust = (PlaneObj(SI4F,SI4B)) ;
                           
                            %                             SelRobust = (PlaneW9(SelIndOffsetF,SelIndOffsetB)).*PlaneW6(SelIndOffsetF,SelIndOffsetB);
                        else
                            SelAlpha = 0.5 ;
                            SelTAlpha=0.5 ; 
                            SelRobust = 0 ;
                        end
                        
                        if sum(isnan(PlaneAlpha(:)))>1
                            a=1 ;
                        elseif sum(isnan(PlaneObj(:)))>1
                            a=1 ;
                        end
                        if sum(SelRobust(:)<0)>1
                            a=1 ;
                        end
                        % generate sample F and B for pixels
                        %                         NumGS =10 ;
                        %                         NumGS4F =round((rand(1)*(NumGS-1)))+1;
                        %                         NumGS4B =round((rand(1)*(NumGS-1)))+1;
                        %
                        %                         TGenSampleF= random(SelModelF,NumGS) ;
                        %                         TGenSampleB= random(SelModelB,NumGS) ;
                        
                        %                         GenSampleF= TGenSampleF (NumGS4F,:);
                        %                         GenSampleB= TGenSampleB (NumGS4B,:);
                        
                       
                        
                        TBlkValNewF(ki,kj,:) = SelF ;
                        TBlkValNewTF(ki,kj,:) = SelTF ;
                        
                        TBlkValNewB(ki,kj,:) = SelB ;
                        TBlkValNewTB(ki,kj,:) = SelTB ;
                        
                        
                        TBlkAlpha(ki,kj) = SelAlpha ; SelAlpha=[] ;
                        TBlkTAlpha(ki,kj) = SelTAlpha ; SelTAlpha=[] ;
                        TBlkRobust(ki,kj) = SelRobust ; SelRobust=[] ;
%%%                                              [ki,kj]
                        % % %                         if (ki==9) && (kj==9)
                        % % %                             a=1 ;
                        % % %                         end
                        % % %
                        % % %                         size(PlaneW1)
                        % % %                         [SelIndOffsetF,SelIndOffsetB]
                        % % %                         a=1
                        TBlkSelWeight(ki,kj,1) = PlaneW1(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,2) = PlaneW2(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,3) = PlaneW3(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,4) = PlaneW4(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,5) = PlaneW5(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,6) = PlaneTW1(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,7) = PlaneTW2(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,8) = PlaneTW3(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,9) = PlaneTW5(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,10) = PlaneCTW1(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,11) = PlaneCTW2(SelIndOffsetF,SelIndOffsetB) ;
                        TBlkSelWeight(ki,kj,12) = PlaneCTW3(SelIndOffsetF,SelIndOffsetB) ;
                        
                       
                        
                        
                    elseif  (TBlkAct(ki,kj)==1)
                        TBlkValNewF(ki,kj,:) = TBlkVal(ki,kj,:);
                        TBlkValNewB(ki,kj,:) = TBlkVal(ki,kj,:);
                        TBlkValNewTF(ki,kj,:) = TBlkTVal(ki,kj,:);
                        TBlkValNewTB(ki,kj,:) = TBlkTVal(ki,kj,:);
                                               
                        TBlkAlpha(ki,kj) = 1 ;
                        TBlkTAlpha(ki,kj)=1 ; 
                        TBlkRobust(ki,kj) = 1 ;
                    elseif (TBlkAct(ki,kj)==5)% pixel is labeled F or labeled B
                        
                        TBlkValNewF(ki,kj,:) = TBlkVal(ki,kj,:);
                        TBlkValNewB(ki,kj,:) = TBlkVal(ki,kj,:);
                        TBlkValNewTF(ki,kj,:) = TBlkTVal(ki,kj,:);
                        TBlkValNewTB(ki,kj,:) = TBlkTVal(ki,kj,:);
                        
                        TBlkAlpha(ki,kj) = 0 ;
                              TBlkTAlpha(ki,kj)=0 ; 
                              
                        TBlkRobust(ki,kj) = 1 ;
                        
                    elseif TBlkAct(ki,kj)==0
                        TBlkValNewF(ki,kj,:) = 255 ;
                        TBlkValNewB(ki,kj,:) = 255 ;
                        
                        TBlkValNewTF(ki,kj,:) = 255 ;
                        TBlkValNewTB(ki,kj,:) = 255 ;
                        
                        TBlkAlpha(ki,kj) = -1 ;
                              TBlkTAlpha(ki,kj)=-1 ; 
                        TBlkRobust(ki,kj) = -2 ;
                        
                    end % EOF TBlkAct(ki,kj)==3
                end  % eof Kj
            end  % eof Ki
            
            IGenF(TLCi : BRCi , TLCj : BRCj,:)=TBlkValNewF ;
            IGenTF(TLCi : BRCi , TLCj : BRCj,:)=TBlkValNewTF ;
            
            
            IGenB(TLCi : BRCi , TLCj : BRCj,:)=TBlkValNewB ;
            IGenTB(TLCi : BRCi , TLCj : BRCj,:)=TBlkValNewTB ;

            
            IRTAlpha(TLCi : BRCi , TLCj : BRCj)=TBlkTAlpha ;
            
            
            IRAlpha(TLCi : BRCi , TLCj : BRCj)=TBlkAlpha ;
            IRobust(TLCi : BRCi , TLCj : BRCj)=TBlkRobust ;
            
            TRIWight(TLCi : BRCi , TLCj : BRCj,:)=TBlkSelWeight ;
            TRIEct(TLCi : BRCi , TLCj : BRCj,1)=E_c ;
            TRIEct(TLCi : BRCi , TLCj : BRCj,2)=E_t ;
            
            
            
            %% =======================================================================================================
            
        elseif (BlkMaskAct(i,j)==-1) % Definit background block
            
            TLCi = (i-1)*BlkH+1 ; TLCj = (j-1)*BlkW+1 ; % block index (Top left corner index) ;
            BRCi = i*BlkH       ; BRCj = j*BlkW   ;     % block index (Bottom right corner) ;
            IGenF(TLCi : BRCi , TLCj : BRCj,:)= 255 ;
            IGenTF(TLCi : BRCi , TLCj : BRCj,:)= 255 ;
            IGenB(TLCi : BRCi , TLCj : BRCj,:)= 255 ;
            IGenTB(TLCi : BRCi , TLCj : BRCj,:)= 255 ;
            
        end % EOF ((BActB(i,j)==0)&& (BActF(i,j)==0))
        
        % get global sample set for interested block ----------------------
        
        
   
        
    end
strcat(num2str((j*100/NumBw)),'% of matting process is completed')
end



%% Return Result ----------------------------------------------------------


RF =IGenF ; % generated Foreground samples
RTF =IGenTF ; % generated Foreground samples
RB =IGenB ; % generated background samples
RTB =IGenTB ; % generated background samples

IRAlpha(MaskAct==1)=1 ; IRAlpha(MaskAct==5)=0 ;
RIAlpha  = IRAlpha ;
RITAlpha  = IRTAlpha ;  % Alpha before optimization
RIRobust = IRobust ;

RIWight = TRIWight  ; % return the effect of weights ;
RIEct= TRIEct ;
%--------------------------------------------------------------------------

