%% in this experiment we test alpha matting based on combined color and
%% texture features .


function  R= WeightedColorTextureMatting(IPath, TrimapPath)

        
        
        IOrg = imread(IPath) ; 
        [IOrgH , IOrgW,Cmode]= size(IOrg) ;
        % Read images and mask act -----------------
%           BlKNumSet = [120 160] ;
          
          BlockSize =[5,5];
          ncol = IOrgW + (mod(IOrgW,BlockSize(2))>0)*(BlockSize(2) - mod(IOrgW,BlockSize(2))) ; 
          nrow = IOrgH + (mod(IOrgH,BlockSize(1))>0)*(BlockSize(1) - mod(IOrgH,BlockSize(1))) ; 
          
         BlKNumSet = [nrow/BlockSize(1),ncol/BlockSize(2)] ;
    
        IOrg = imresize(IOrg,[nrow ncol]); IOrg=uint8(IOrg);
        TSize =1;
        
        I=IOrg  ; I= imresize(I,TSize);I = double(I) ;
        [Ih,Iw,Cmod]= size(I) ;
        IGray=rgb2gray(uint8(I));IGray=double(IGray) ;  % Gray Level Image
        
        IGrad = abs(gradient(double(IGray))) ;
        IGradRGB (:,:,1) = abs(gradient(double(I(:,:,1)))) ;
        IGradRGB (:,:,2) = abs(gradient(double(I(:,:,2)))) ;
        IGradRGB (:,:,3) = abs(gradient(double(I(:,:,3)))) ;
        
        CutThr = .04 ;

        % Trimap and MaskAct ----------
        ITrimap= imread(TrimapPath);
        ITrimap=ITrimap(:,:,1) ; 
        ITrimap = imresize(ITrimap,[nrow ncol]); ITrimap=uint8(ITrimap);
        ITrimap= imresize(ITrimap,TSize);
        RmaskB= (ITrimap < 50);  RmaskF= (ITrimap >200 );
        % Extrapolate maskes ----
        RMDilaterFilter = ones(3)>0 ;
        RmaskB = imdilate(RmaskB,RMDilaterFilter ) ;
        RmaskF = imdilate(RmaskF,RMDilaterFilter ) ;
        RmaskBOrg=RmaskB ;
        RmaskFOrg=RmaskF ;
        MaskActOrg = RmaskBOrg * 5 + RmaskFOrg ;   % get the marked F&B pixels
        MaskActOrg(MaskActOrg==0)=3 ;
        MaskAct=MaskActOrg ;
        MaskAct(MaskAct==0)=3 ;
        % LabelExpansion -------------------------------------------------
        %         ExpDist = 7 ;
        %         ExpThr=8/256 ;
        %         TCloseStrl= ones(5) ;
        %         [RMaskFExp ,RMaskBExp] = LabelExpansion (I, RmaskF,RmaskB , ExpDist, ExpThr)  ;
        %         RmaskF = imclose(RMaskFExp,TCloseStrl) ;
        %         RmaskB = imclose(RMaskBExp,TCloseStrl) ;
        %         MaskAct = zeros (Ih,Iw) ; % 0= inactive , 1= definit foreground ,% 5= definit background , 3= potential pixels
        %         MaskAct = RmaskB * 5 + RmaskF ;   % get the marked F&B pixels
        
        % Color Format of Image (RGB or Ycbcr ) ---------------------------
        ColFormat= 'rgb' ;
        %ColFormat= 'ycbcr' ;
        if strcmpi(ColFormat,'ycbcr')
            IColor=double(rgb2ycbcr(uint8(I))) ;
        elseif strcmpi(ColFormat,'rgb')
            IColor = double(I) ;
        end
        
        
        
        

        
           

            
            MaskActWOExp = abs(MaskActOrg -MaskAct) ; MaskActWOExp(MaskActWOExp>0)= 3 ;
            MaskActWOExp(RmaskFOrg)=1 ; MaskActWOExp(RmaskBOrg)=5 ;
            
            
            %% 1 ) Alpha Matting --------------------------------------------
            
            % 1-1) Texture Image Construction =============================
            TextMode = 2 ; % 1 for Gray , , 2 for Corlo Texture
            TErodLen = 200 ; % length of borders pixels
            
            TRSize=3 ;
            IColor2=imresize(IColor, TRSize) ;
            IGray2= double(rgb2gray(uint8(IColor2)));
            MaskAct2 =zeros(size(IGray2))+3 ;
            
            MaskAct2( imresize(RmaskF,TRSize)) =1 ;
            MaskAct2( imresize(RmaskB,TRSize)) =5 ;
            
            NumWLevel=2;
            TSize = 2^(NumWLevel-1) ;
            %             %             tic
            [TextureImage  , TAlphaCF]  = Fun_TextureImageConstructionV18 (IColor2 ,IGray2, MaskAct2 ,TErodLen,TextMode, NumWLevel ); % sets are not  equalized .
            toc
            ITexture = double(imresize(TextureImage,1/TRSize )) ;
            TextureImage= ITexture ;
            IGradTex = abs(gradient(double(rgb2gray(uint8(ITexture))))) ;
            IGradTex =  FunCoefHistCut_V2 (IGradTex,CutThr ) ;
            
            % Texture Effect of Mask Expansion  ---------------------------
            ExpDist = 10 ;     ExpThr=3/255 ;
            ExpDist2= 10 ;     ExpThr2 = 3/255 ;
            TCloseStrl= ones(3) ;
            [RMaskFExp ,RMaskBExp] = LabelExpansion (I, RmaskF,RmaskB , ExpDist, ExpThr)  ;
            [RMaskFExpTex ,RMaskBExpTex] = LabelExpansion (ITexture, RmaskF,RmaskB , ExpDist2, ExpThr2)  ;
            
            RMaskFExp2 = (RMaskFExp|RMaskFExpTex) ;
            RMaskBExp2 = (RMaskBExp|RMaskBExpTex) ;
            RmaskF = imclose(RMaskFExp2,TCloseStrl) ;
            RmaskB = imclose(RMaskBExp2,TCloseStrl) ;
            MaskAct = zeros (Ih,Iw) ; % 0= inactive , 1= definit foreground ,% 5= definit background , 3= potential pixels
            MaskAct = RmaskB * 5 + RmaskF ;
            MaskAct(MaskAct==0)=3 ;
            
            ITrimap=zeros(size(MaskAct))+128 ;
            ITrimap(MaskAct==1)=255 ;
            ITrimap(MaskAct==5)=0 ;
            
            
            
            % 1-2) Color and Texture Matting --------------------------------
          
%               BlKNumSet = [60 80] ;
            EWC= [ 0 1 0 0 0 ] ;  % effectiveness weights for color in Matting Process
            %           EWC= [.5 1 1 0 .5 ] ;
            EWT= [0 0 0 0 ] ;  % effectiveness weights for texture in Matting Process
            EWCT= [ 0.5 0 0 ] ;  % effectiveness weights for color and texture conjunction  in Matting Process
            NumKg= 4 ;
            NumS4C = 2; % number of samples generated for each clusster inside the blocks
            
            msg = ' matting  process is started.'             
            [RAlpha,RConf, RTexAlpha ,RIFG,RIBG,RIWeight, RIEct ] = Fun_SharedDistMatting_V55_7_3 (IColor,ITexture, IGrad , MaskAct,  BlKNumSet,NumKg, EWC,EWT,EWCT,NumS4C) ;
            msg = ' matting process is completed.' 
            msg = ' refinement process is started.' 
            % 1-3 ) Refinement Step using laplacianmatrix ==============
            RC= RIWeight(:,:,2).*RIWeight(:,:,10) ;RConf=RC ;
            RConf(RmaskF)=1 ; RConf(RmaskB)=1 ;
            RConf=(RConf) ;
            pack=[] ;
            pack(:,:,1) = uint8(RAlpha*255 ) ;
            pack(:,:,2) = uint8((RConf)*255 ) ;
            pack(:,:,3) = uint8(ITrimap) ;
            alpha = Fun_SmoothingMatting(IColor, pack) ;
            msg = ' refinement process is completed.' 
            
            
            %% 3) Save Results  -----------------------------------------------
            
            NewAlpha= alpha ; NewAlpha(NewAlpha<0)=0 ; NewAlpha(NewAlpha>1)=1 ; NewAlpha=(NewAlpha*255) ;
            NewAlpha= uint8(imresize(NewAlpha,[IOrgH,IOrgW])) ;
            imwrite(NewAlpha,'EstimatedAlpha.png' ) ;  % save color estimated matte   
            
            R= NewAlpha ; 
           
         %%ReturnResult    
         
%          RefinedAlpha = alpha ; 
%          TextureImage = TextureImage ; 
         
            
            
end  %% End Function :) ===================================================
% =========================================================================