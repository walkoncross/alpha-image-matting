function [RC,RLBL,RCriteria] = Fun_Find_NumCluster(Data,MaxK, mode)
%% in this function number of optimum cluster is computed base of simple
%% criteria like elbow and AIC or BIC
%
%--------------------------------------------------------------------------
%%  Input Variables :
%       Data is NxP data matrix where P is number of featurs
%       MaxK    : indicates the maximum nuber of clussters that should be evaluate
%       mode  0 : use elbow criterion to find optimal number of clussters
%             1 : use AIC criterion to find optimal number of clussters
%             2 :
%
%
%% ------------------------------------------------------------------------



if mode ==0  % use elbow criterion and Kmean to find optimal number of clussters
    
    %     MaxK = 10 ;
    Criteria= zeros(1,MaxK); % Criteria of clusters
    IDXSet = zeros(size(Data,1),MaxK);
    
    for k=2 : MaxK
        
        
        [IDX,C,sumd,D] = kmeans(Data,k,'emptyaction','singleton') ;
        NumSamp=[] ;
        for i=1 : k
            NumSamp(i)=sum(IDX==i) ;
        end
        NumSamp(NumSamp==0)= 1 ;
        SWithin =  sum(sumd'./NumSamp) ;  % variance withing clusters
        
        a= dist(C');
        SBetween = sum(a(:))/2  ; % variance between clusters
        
        Criteria1(k)=(SWithin/SBetween);
        Criteria(k)= (k^2)*(SWithin/SBetween);
        
        IDXSet(:,k)=IDX ;
        %    Criteria1(k)= (SWithin/SBetween)
        
    end
    
    %
    % [TminK,Kind]=min(Criteria) ;
    
    a=Criteria1(2:end)- [Criteria1(3:end),Criteria1(end)]
    a= [100 , a] ;
    
    [TmaxK,Kind]=max(a(a<=1.5)) ;
    
    if length(TmaxK)==0
        [TmaxK,Kind]=min(Criteria1) ;
    end
    CriteriaInd = (1:1:MaxK);
    % Kind = min(CriteriaInd(Criteria1<=TmaxK));
    Kind = min(CriteriaInd(a<=TmaxK));
    
    %% Return Result ----------------------------------------------------------
    RC =Kind ; %Return number of optimal clusters
    RLBL = IDXSet(:,Kind) ;
    RCriteria= Criteria1 ;
    
elseif mode==1 % AIC criterion to find optimal number of clussters using Gaussian Mixture Models
    
    
    %   MaxK = 10 ;
    
    AIC_Criteria = zeros(1,MaxK);
    obj = cell(1,MaxK);
    for k=2 : MaxK
        
        obj{k} = gmdistribution.fit(Data,k,'Regularize',0.0005);
        AIC_Criteria(k)= obj{k}.AIC;
        
        k
        
    end
    
    
    
    
    
    
    
    % preprocessing to improve the AIC criteria ---------------------------
    a=AIC_Criteria(2) ;
    
    for k=2 :MaxK
        if AIC_Criteria(k)<a
            
            a=AIC_Criteria(k)
        else
            AIC_Criteria(k)= a ;
        end
        
        
    end
    %----------------------------------------------------------------------
    % define a way to find optimal number of clussters --------------------
        % use chart TO FIND BEST POINT
    %     |              + + + + +
    %     |            +
    %     |          +
    %     |        +
    %     |      +
    %     |     +
    %     |     +
    %     |    +
    %     |  +
    %     |--------------------------->
    a= AIC_Criteria ;
    
    AIC_dist = AIC_Criteria (1: MaxK-1)- AIC_Criteria(2:MaxK);
    AIC_dist=[0 , AIC_dist] ;
    AcumDist = cumsum(AIC_dist) ; AcumDist =AcumDist/AcumDist(MaxK);
    
    chart(1,:)=AcumDist ;
    chart(2,:)= [1:MaxK]/MaxK ;
    

    GoalPoint=[1;0] ;
    TempChartEval= sqrt(10*(chart(1,:)-GoalPoint(1)).^2 + (chart(2,:)-GoalPoint(2)).^2);
    
    [minVal , OPNumC] = min(TempChartEval) ; %OPNumC is optimal number of clusters
    
    % get clusster index of data ------------------------------------------
    DataClusterIndex = cluster(obj{OPNumC},Data) ;
   
        %% Return Result ----------------------------------------------------------
    RC =OPNumC ; %Return number of optimal clusters
    RLBL = DataClusterIndex ; % label of data points 
    RCriteria= TempChartEval ;
    
    
end % endof mode

