function RSmoothAlpha = ComprehensiveSamplingMatting (IOrg, ITrimap)

        

        [IOrgH , IOrgW,Cmode]= size(IOrg) ;
        % Read images and mask act -----------------
        BlkHlen = 1 ;  BlkWlen = 1 ;
        BlkSize=[BlkHlen ,BlkWlen ] ;
        nrow =IOrgH + (BlkHlen- (IOrgH-floor(IOrgH / BlkHlen)*BlkHlen))*((IOrgH-floor(IOrgH / BlkHlen)*BlkHlen)>0);
        ncol =IOrgW + (BlkWlen- (IOrgW-floor(IOrgW / BlkWlen)*BlkWlen))*((IOrgW-floor(IOrgW / BlkWlen)*BlkWlen)>0) ;
        
        
               IOrg = imresize(IOrg,[nrow ncol]); IOrg=uint8(IOrg);
        TSize =1;

        I=IOrg  ; I= imresize(I,TSize);I = double(I) ;
        [Ih,Iw,Cmod]= size(I) ;
        IGray=rgb2gray(uint8(I));IGray=double(IGray) ;  % Gray Level Image


        % Trimap and MaskAct ----------
        ITrimap = imresize(ITrimap,[nrow ncol]); ITrimap=uint8(ITrimap);
        ITrimap= imresize(ITrimap,TSize);
        RmaskB= (ITrimap < 50);  RmaskF= (ITrimap >200 );
        % Extrapolate maskes ----
        RMDilaterFilter = ones(3)>0 ;
%         RmaskB = imdilate(RmaskB,RMDilaterFilter ) ;
%         RmaskF = imdilate(RmaskF,RMDilaterFilter ) ;
        RmaskBOrg=RmaskB ;
        RmaskFOrg=RmaskF ;
        MaskActOrg = RmaskBOrg * 5 + RmaskFOrg ;   % get the marked F&B pixels
        MaskActOrg(MaskActOrg==0)=3 ;
        MaskAct=MaskActOrg ;
        MaskAct(MaskAct==0)=3 ;
        % LabelExpansion -------------------------------------------------
%         ExpDist = 7 ;
%         ExpThr=8/256 ;
        RMaskFExp=RmaskF ; RMaskBExp=RmaskB ;
        ExpThr_U=9/256 ;
        ExpThr_D=1/256 ;
        ExpThrDist = ExpThr_U-ExpThr_D ;
        MaxIterExp=9 ;
        'Trimap expansion ...'
        for i=1 : MaxIterExp            
            ExpDist =i ;
            ExpThr =ExpThr_U - i* ExpThrDist / MaxIterExp ;
            [RMaskFExp ,RMaskBExp] = LabelExpansion (I, RMaskFExp,RMaskBExp , ExpDist, ExpThr)  ;
            
        end
%         TCloseStrl= ones(2) ;
%         RMaskF = imerode(RMaskFExp,TCloseStrl) ;
%         RMaskB = imerode(RMaskBExp,TCloseStrl) ;
       RmaskF =  RMaskFExp ; 
       RmaskB =  RMaskBExp ; 

       MaskAct = zeros (Ih,Iw) ; % 0= inactive , 1= definit foreground ,% 5= definit background , 3= potential pixels
       MaskAct = RmaskB * 5 + RmaskF ;   % get the marked F&B pixels
       MaskAct(MaskAct==0)=3 ;


            % -------------------------------------------------------------
            MaxLevel=6 ; 
            BThickness=[5 15] ; 
            for l=1 : MaxLevel 
                if l>2 
            BThickness(l) =2* BThickness(l-1) +BThickness(l-2)  ;
                   
                end
%             RBoundaryF =RmaskF - imerode(RmaskF, ones(BThickness)>0 );
%             RBoundaryB =RmaskB - imerode(RmaskB, ones(BThickness)>0 );
%             MaskActBoundary= 5*RBoundaryB + RBoundaryF ; MaskActBoundary(MaskAct==3)=3 ;    
IColor=I ; 
            LSmoothness = [3 3 5] ;
            
            mode = 'CS' ; CutThr = 0.0001 ;
            tic
            [L14F,L14B ,L24F, L24B  ]= Fun_MultiLevelPixelBoundarySmpSel_Color_V1(IColor,  MaskAct,  BThickness(l),LSmoothness,mode,CutThr) ;
            toc
            ActL4F{1,l} = L24F ;
            ActL4B{1,l} = L24B ;
            
            end
            
           %**************************************************************
           % UnKnown Pixel Level Correspondency 
            ULbl2F = zeros(Ih,Iw)+MaxLevel ; 
            ULbl2B = zeros(Ih,Iw)+MaxLevel ; 
            
           for l= MaxLevel :-1:1
           TRmaskF1 = imdilate(RmaskF, ones(BThickness(l))>0 );  
           TRmaskB1 = imdilate(RmaskB, ones(BThickness(l))>0 );  
           
           ULbl2F(TRmaskF1)= l ; 
           ULbl2B(TRmaskB1)= l ; 
               
           end
             ULbl2F(RmaskF)=0 ; ULbl2F(RmaskB)=0 ; 
             ULbl2B(RmaskF)=0 ; ULbl2B(RmaskB)=0 ; 
           %***************************************************************
             
             
             
             
            % *************************************************************
            % 1-2) Color and Texture Matting ------------------------------
            EWC = [0.5 0.5 0  2  ] ;
            SelMod=1 ;   % 1) mean values for image .    2)  random samples per pixel
            Normal=0 ;
            tic
            SmpMod= 'mean' ;
            [TR_FV_F , TR_FV_B,TR_FVAlpha , TRobust,TRSelWeight,TRICohenD]=Best_FBSelection_MLevelClusstering (IColor, MaskAct,  ActL4F, ActL4B, ULbl2F, ULbl2B , EWC , SmpMod);
            
            RAlpha= TR_FVAlpha ; 
            % 1-3 ) Refinement Step using laplacianmatrix ==============
            RConf1=sqrt(TRobust) ; 
            RConf1(RmaskF)=1 ; RConf1(RmaskB)=1 ;
            RConf1=(RConf1) ;
            pack=[] ;
            pack(:,:,1) = uint8(TR_FVAlpha*255 ) ;
            pack(:,:,2) = uint8((RConf1)*255 ) ;
            pack(:,:,3) = uint8(ITrimap) ;
            alpha = Fun_SmoothingMatting(IColor, pack) ;
            
            
   
            
          %% Return The results -------------------------------------------
          
               RSmoothAlpha = alpha ; 

          
          % ---------------------------------------------------------------
            
            
    end