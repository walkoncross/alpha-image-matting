function R= FunCoefHistCut_V4 (Coef, CutThr)
%% the obly difference with version2 is that here We return Original Values Here -------------------------------------
% in this function we cut the head and tail of histogram and expand it
[Coef_H , Coef_W , Cmod]= size(Coef) ; 

if Cmod==1 

Coef= Coef - min(Coef(:)) ; 
a= min(Coef(:)) ;  b=max(Coef(:)) ; 
Coef=Coef.*255/ max(Coef(:)) ;
CoefHist = imhist(uint8(Coef)) ; 
CoefHistR = (CoefHist(end:-1:1)) ; % reversed histogram

HistCum= cumsum(CoefHist);
HistCumR= cumsum(CoefHistR);  % commulative sum of reversed Histogram

HistCheck= (HistCum > HistCum(end)*CutThr); 
HistCheckR= (HistCumR > HistCum(end)*CutThr); 

[TVal CBottomThr]=max(HistCheck) ;  % find threshold for low color in histogram
[TVal CTopThr]=max(HistCheckR) ; 
CTopThr = 255  - CTopThr ;  % find threshold for top colors values


Coef(Coef >= CTopThr)=CTopThr; Coef(Coef <= CBottomThr)=CBottomThr;
Coef= Coef - min(Coef(:)) ; Coef=Coef.*255/ max(Coef(:)) ; 
% bar(imhist(uint8(Coef)))--
% Return Result ------------------------------
% Return To original Value
Coef= Coef * (b/255) ; Coef = Coef+a ;  
R = Coef ;

elseif Cmod==3 
    for i=1 : Cmod
     R(:,:,i)= FunCoefHistCut_V4 (Coef(:,:,i), CutThr) ;
     
    end
    
else
    
    error ('Image should be 1 dimension or 3 dimensions')
    
end

