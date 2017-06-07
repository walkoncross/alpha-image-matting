function R = Fun_ComprehensiveMatting(IOrg, ITrimap) 


        
       [IOrgH , IOrgW,Cmode]= size(IOrg) ;
        % Read images and mask act -----------------
        BlkHlen = 1 ;  BlkWlen = 1 ;
        BlkSize=[BlkHlen ,BlkWlen ] ;
        nrow =IOrgH + (BlkHlen- (IOrgH-floor(IOrgH / BlkHlen)*BlkHlen))*((IOrgH-floor(IOrgH / BlkHlen)*BlkHlen)>0);
        ncol =IOrgW + (BlkWlen- (IOrgW-floor(IOrgW / BlkWlen)*BlkWlen))*((IOrgW-floor(IOrgW / BlkWlen)*BlkWlen)>0) ;
        
        
        IOrg= imread(filename1) ; [IOrgH , IOrgW,Cmode]= size(IOrg) ;
        IOrg = imresize(IOrg,[nrow ncol]); IOrg=uint8(IOrg);
        TSize =1;
        %         for d=1 : 3
        %             IOrg(:,:,d)= histeq(uint8(IOrg(:,:,d)));
        %         end
        I=IOrg  ; I= imresize(I,TSize);I = double(I) ;
        [Ih,Iw,Cmod]= size(I) ;
        IGray=rgb2gray(uint8(I));IGray=double(IGray) ;  % Gray Level Image        
                
       
        CutThr = .04 ;
        

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
        for i=1 : MaxIterExp            
            ExpDist =i ;
            ExpThr =ExpThr_U - i* ExpThrDist / MaxIterExp ;
            [RMaskFExp ,RMaskBExp] = LabelExpansion (I, RMaskFExp,RMaskBExp , ExpDist, ExpThr)  ;
            i
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
            MaxLevel= 9  ;
            BThickness=[5 15] ;
            for l=1 : MaxLevel
                if l>2
                    BThickness(l) =BThickness(l-1) +BThickness(l-2)  ;
                    
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
            [TR_FV_F , TR_FV_B,TR_FVAlpha , TRobust,TRSelWeight,TRICohenD,L2TR_FVAlpha , L2TRobust]=Best_FBSelection_MLevelClusstering_V1_1 (IColor, MaskAct,  ActL4F, ActL4B, ULbl2F, ULbl2B , EWC , SmpMod);
            
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
            


            %% 3) Save Results  -----------------------------------------------
            IterInd=1  ;
            WRoot ='.\' ;
            
            WData = strcat(WRoot,FnameP1(1:end-4),'_',TrimapSet{TrInd},'_Iter',num2str(IterInd),'_Data.mat') ;
            WRC = strcat(WRoot,TrimapSet{TrInd},'\','Raw_RA_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRT = strcat(WRoot,TrimapSet{TrInd},'\','Raw_RTA_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRCol = strcat(WRoot,TrimapSet{TrInd},'\','L0_RLapAlpa_Color_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRFirstSmoothAlpha = strcat(WRoot,TrimapSet{TrInd},'\','FirstSmoothAlpha_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRFirstSmoothAlpha2 = strcat(WRoot,TrimapSet{TrInd},'\','L2FirstSmoothAlpha_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRAlphaC = strcat(WRoot,TrimapSet{TrInd},'\','L1_RLapAlpa_Color_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRAlphaT = strcat(WRoot,TrimapSet{TrInd},'\','L1_RLapAlpa_Texture_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRIE_c = strcat(WRoot,TrimapSet{TrInd},'\','IEc_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRIE_t = strcat(WRoot,TrimapSet{TrInd},'\','IEt_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRAlphaFinalOldV = strcat(WRoot,TrimapSet{TrInd},'\','L2_RLapAlpa_FinalOldV_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRCol2 = strcat(WRoot,TrimapSet{TrInd},'\','L1_RSmoothFunAlpa_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRCol3 = strcat(WRoot,TrimapSet{TrInd},'\','L1_RSmoothFunAlpaOldV_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRTImg = strcat(WRoot,TrimapSet{TrInd},'\','TextureImage_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRTMaskAct = strcat(WRoot,TrimapSet{TrInd},'\','MaskAct_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRTConf = strcat(WRoot,TrimapSet{TrInd},'\','L0RConfidence_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            WRTConf2 = strcat(WRoot,TrimapSet{TrInd},'\','L2RConfidence_',FnameP1(1:end-4),'_Iter',num2str(IterInd),'.png') ;
            
            NewAlpha= alpha ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRFirstSmoothAlpha ) ;  % save color estimated matte
            
            
            NewAlpha= alpha2 ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRFirstSmoothAlpha2 ) ;  % save color estimated matte
            
            
            
            
            % Save Combined Alpha -----------------------------------------
%             NewAlpha= RIEct(:,:,1) ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
%             NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
%             imwrite(NewAlpha,WRIE_c ) ;  % save color estimated matte
%             
            % Save Combined Alpha -----------------------------------------
%             NewAlpha= RIEct(:,:,2) ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
%             NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
%             imwrite(NewAlpha,WRIE_t ) ;  % save color estimated matte
            
            % Save Combined Alpha -----------------------------------------
            NewAlpha= RAlpha ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRC ) ;  % save color estimated matte
            
            

            % Save Just Color Alpha ---------------------------------------
            NewAlpha= alpha ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRCol ) ;  % save texture estimated matte
            
            NewAlpha= RConf ;   NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRTConf2 ) ;  % save confidence values
            
            
            NewAlpha= RConf1 ;   NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRTConf ) ;  % save confidence values
            

            
            

            
            
            NewAlpha= MaskAct ;   NewAlpha=(NewAlpha*51) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,WRTMaskAct ) ;  % save confidence values
            
            
%             save(WData) % save Data
            
            % -------------------------------------------------------------
          
          
            
            
















end