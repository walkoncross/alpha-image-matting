function R = GetNextLevelCluster(L0 , Nfeatures , ClusterMod,LSmoothness)


Cmod = Nfeatures(1) ;
Tmod = Nfeatures(2) ;
FVmod= Cmod+Tmod; 





L0NumC = size(L0,1)  ;
InsertCluster=0 ;
for ci=1 : L0NumC
    
    
    L0TSmp_mean = L0{ci,1};
    L0TSmp_std = L0{ci,2};
    L0TSmp = L0{ci,3};
    L0TNumSmp = L0{ci,4};
    L0TSmp_2DInd = L0{ci,5};
    L0TSmp_LInd = [1:1 : L0TNumSmp]';
    L0TSmp_GInd = L0{ci,6};
    if L0TNumSmp>0
    if strcmpi(ClusterMod,'color')
        
        for chi=1 : Cmod
            [TNumP(chi,1),THist(chi,:)]= FunGetBlockPeak (L0TSmp(:,chi) ,ones(L0TNumSmp,1)>0, LSmoothness );
        end        
        L1NumC = round(max(TNumP)); % estimated number of clussters ----------
           L1NumC = min(round(L0TNumSmp/2),L1NumC) ;
        
        TKmeanFV= L0TSmp(:,1:Cmod) ; 
        
    
    ci
    
    elseif strcmpi(ClusterMod,'texture')
        
        for chi=1 : Tmod
            [TNumP(chi,1),THist(chi,:)]= FunGetBlockPeak (L0TSmp(:,chi+Cmod) ,ones(L0TNumSmp,1)>0, LSmoothness );
        end
        L1NumC = round(max(TNumP)); % estimated number of clussters -------
           L1NumC = min(round(L0TNumSmp/2),L1NumC) ;
        TKmeanFV= L0TSmp(:,Cmod+1:end) ;
        
    elseif strcmpi(ClusterMod,'spatial')
        TempInd= L0TSmp_2DInd; 
        TempInd =TempInd- repmat(min(TempInd),[L0TNumSmp,1]); 
        TempInd =double(uint8(TempInd*255./repmat(max(TempInd),[L0TNumSmp,1]))); 
       ci
        for chi=1 : 2
            [TNumP(chi,1),THist(chi,:)]= FunGetBlockPeak (TempInd(:,chi) ,ones(L0TNumSmp,1)>0, LSmoothness );
        end
        L1NumC = round(max(TNumP)); % estimated number of clussters -------
        L1NumC = min(round(L0TNumSmp/2),L1NumC) ;
        TKmeanFV= L0TSmp_2DInd+ rand(size(L0TSmp_2DInd))/100 ;
        
    end  % EOF IF Color Mod

        if ((L0TNumSmp*.05) > L1NumC )
            
            obj = gmdistribution.fit(TKmeanFV,L1NumC,'Regularize',eps,'SharedCov',true , 'CovType','diagonal');
            [L1_Idx,nlogl] = cluster(obj,TKmeanFV) ;
%             [L1_Idx , L1_Ctrs] = kmeans(TKmeanFV,L1NumC,'start','cluster','emptyaction','singleton');
        else
%             [L1_Idx , L1_Ctrs] = kmeans(TKmeanFV,L1NumC,'emptyaction','singleton');
              obj = gmdistribution.fit(TKmeanFV,L1NumC,'Regularize',eps,'SharedCov',true , 'CovType','diagonal');
            [L1_Idx,nlogl] = cluster(obj,TKmeanFV) ;
            
        end

   
    
    
    for L1_ci=1 : L1NumC
        
        L1TSmp=[] ; TSmplGlobalInd1D=[] ; TSmpl2DInd=[] ;
        
        L1TLInd  = L0TSmp_LInd(L1_Idx==L1_ci) ; % get local index of clusster members
        
        for chi=1 : FVmod
            L1TSmp(:,chi)=L0TSmp(L1TLInd,chi);
            
        end
        
        NumL1TSamp = size(L1TLInd,1) ;
        TCSmpCov = cov(L1TSmp) ;
        
        TCSmpStd = std(L1TSmp) ;
        TCSmpMean= mean(L1TSmp) ;
        TCPrior = NumL1TSamp/L0TNumSmp ;
        
        
        TSmpl2DInd(:,1)=L0TSmp_2DInd(L1TLInd,1);
        TSmpl2DInd(:,2)=L0TSmp_2DInd(L1TLInd,2);
        
        % TSmpl1DInd=FSet_Ind1D(TCInd);
        
        TSmplGlobalInd1D = L0TSmp_GInd(L1TLInd) ;
        
        % -----------------------------------------------------------------
        InsertCluster=InsertCluster+1 ; 
        L1{InsertCluster,1}=TCSmpMean ;  % 1D Index of samples
        L1{InsertCluster,2}=TCSmpStd ;  % 1D Index of samples
        
        L1{InsertCluster,3}=L1TSmp ;  % Color Features of samples
        L1{InsertCluster,4}=NumL1TSamp;  % Color Features of samples
        
        L1{InsertCluster,5}=TSmpl2DInd ;    % 2D Index of samples
        L1{InsertCluster,6}=TSmplGlobalInd1D ;  % 1D Index of samples
        L1{InsertCluster,7}=TCSmpCov ;  % Covariance of samples 
        % -----------------------------------------------------------------
    
        
    end % Eof L1_ci
    
    end
    
        
        
    
end % Eof L0 Ci


%% Return Result 

R= L1 ; % return clussters  :)


end % EOF Function --------------------------------------------------------