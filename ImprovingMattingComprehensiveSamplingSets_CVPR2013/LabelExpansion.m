function [RmaskFExp ,RmaskBExp] = LabelExpansion (I, RmaskF,RmaskB , ExpDist, ExpThr) 
%% this function expand the foreground and background labels to unknown area . 
%% the expansion is based on color similarity  and when unknown pixels are
%% close in color space less than ExpThr , the label is assinged to them
%% based on theirspatial clossness and color similarities





%%






I = double(I)/255;
[Ih, Iw, Cmod] =  size(I) ; 



RmaskR = (RmaskB | RmaskF ) ; 
% RmaskU=(~RmaskR) ;
RmaskR2 = double(RmaskB*5) +double(RmaskF)  ;
RmaskRExpand = double(RmaskB*5) +double(RmaskF)  ; 

% Activation Mask -------------------
% we use this process to reduce computational time ------------------------
DilateFilt= ones(2*ExpDist) ;
RMF_Dilate = imdilate(RmaskF,DilateFilt) ; 
RMB_Dilate = imdilate(RmaskB,DilateFilt) ; 

RmaskU= ((RMF_Dilate - RmaskF) + (RMB_Dilate - RmaskB))>0 ; 
% -------------------------------------------------------------------------



Ki =ExpDist ;
Kc= ExpThr ; 

Ci = (2*Ki)/2 + 1 ; Cj = (2*Ki)/2 +1 ;
for i=1 : 2*Ki+1
    for j=1 : 2*Ki+1
        DistPlane  (i,j)= sqrt((i-Ci)^2 + (j-Cj)^2) ;
        
    end
end

for i=1 : Ih
    for j=1 : Iw
        
        if RmaskU (i,j)==1
            TBlockC = [] ;    TBlockR = [] ;   TColDiff=[] ;  TCDistortion=[] ;
            
            TRCi =max( i- Ki , 1 ) ;  TRCj =max( j- Ki , 1 )  ;
            BLCi =min (i+Ki , Ih)  ; BLCj =min (j+Ki , Iw)   ;
            
            DTRCi =Ci-( i - TRCi) ;  DTRCj =Cj-(j - TRCj)  ;
            DBLCi =Ci+( BLCi-i) ; DBLCj = Cj+(BLCj-j ) ;
            
            PCVal (1) = I(i,j,1); PCVal (2) = I(i,j,2);PCVal (3) = I(i,j,3);
            
            TBlockC = I (TRCi:BLCi , TRCj :BLCj , :) ;
            TBlockR = RmaskR (TRCi:BLCi , TRCj :BLCj , :) ;
            TBlockFBUMask = RmaskR2 (TRCi:BLCi , TRCj :BLCj ) ;
            TBlockDist = DistPlane (DTRCi:DBLCi , DTRCj :DBLCj ) ;
            % check color distortion  ;
            TColDiff(:,:,1) = TBlockC (:,:,1) - PCVal(1) ;
            TColDiff(:,:,2) = TBlockC (:,:,2) - PCVal(2) ;
            TColDiff(:,:,3) = TBlockC (:,:,3) - PCVal(3) ;
            TColDist = sqrt(sum(TColDiff.^2,3)) ;
%             [TColDistMin , TColDistMinInd] = min(TColDist(:)) ;
%             TChk2= (TColDist == TColDistMin) ;
            TCDistortion = TColDist < Kc ;
            
            TChk1= (TCDistortion .* TBlockR)>0 ;
            if  sum (TChk1(:)) >= 1
                TD1= TBlockDist(TChk1) ; 
                [TD1min, TDminInd]=min(TD1) ; 
                TSelPixMask = TBlockFBUMask(TChk1) ;
                
                
                RmaskRExpand (i,j)= TSelPixMask(TDminInd) ;
            end
            
            
        end
        
        
    end
end





%% Return Result ----------------------------------------------------------

RmaskFExp = (RmaskRExpand==1) ;
RmaskBExp = (RmaskRExpand==5) ;



%% ------------------------------------------------------------------------
% 
% subplot(1,2,1) ; imagesc(RmaskR2); figure(gcf); title('Trimap')
% subplot(1,2,2) ; imagesc(RmaskRExpand); figure(gcf); title('Expansion of Known Regions')
% a =(RmaskRExpand~=0)&(RmaskR2==0 ) ;

