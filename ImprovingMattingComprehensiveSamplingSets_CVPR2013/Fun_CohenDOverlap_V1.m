function R = Fun_CohenDOverlap_V1 (MuB, SigmaB, NumSmp4B , MuF , SigmaF,NumSmp4F, Mod ) 
% the V3 is similar to V2 
% In version V3 the Cohend method is used to find intersection between
% distribusions  .
% ------------------------- ------------------------------------------------

[NumB , Cmod ] = size(MuB) ; 
[NumF , Cmod ] = size(MuF) ; 
OverlapWeight = zeros(NumF,NumB) ; %weight of overlapping distributions



%%Cohen D methods 
for i=1 : NumB 
    for j=1 : NumF 
   TMUB = MuB(i,:);     
   TStdB= SigmaB(i,:) ; 
   
   TMUF = MuF(j,:);
   TStdF= SigmaF(j,:) ; 
   
   TMuDist = abs(TMUF - TMUB) ; 
%    TS = (TStdB + TStdF)/2 ; 
%    Td =mean( TMuDist./TS );
   
   TS2 = ((NumSmp4B(i)-1)* (TStdB.^2) + (NumSmp4F(j)-1)* (TStdF.^2)) ;
   TS2 = sqrt(TS2 / (NumSmp4B(i)+NumSmp4F(j)-2));
   Td2 =mean( TMuDist./TS2 );
   
%    OverlapWeight(j,i)= exp(-1/Td);
if  strcmp(Mod,'Org')   
    
% OverlapWeight(j,i)= exp(-1/Td2);  
OverlapWeight(j,i)=1/Td2 ; 
elseif strcmpi(Mod,'CohenD')
OverlapWeight(j,i)=Td2 ;     
    
elseif strcmp(Mod,'exp')
    
    
    
    
OverlapWeight(j,i)= 1-exp(-1/Td2);        

elseif strcmpi(Mod,'CohenD')
    
OverlapWeight(j,i)= Td2;            
end

   

        
        
    end
end
    
   %% Return Result ------------
   if  strcmp(Mod,'log')              
%        OverlapWeight=OverlapWeight./ max(OverlapWeight(:));       
   elseif strcmp(Mod,'exp')
%        OverlapWeight=OverlapWeight./ max(OverlapWeight(:));              
   end

   R= OverlapWeight ;
   
% ---------- --------------------
    