


function [RPeakNumB,RHist]= FunGetBlockPeak (I ,RmaskB,varargin)
% compute the number of peak in the histogram of I ; 

% varargin{1} is smooth size for histogram 

HistSet=[] ; 

[Ih, Iw,Cmod]= size(I) ;

for channeli=1 : Cmod

    IColor = I(:,:,channeli) ;

Bind = 0  ;  % number of F&B marked Pixels
PeakNumB=0 ;
%  1-1) get foreground and background marked pixels--------
% for i=1 : Ih
%     for j=1 : Iw
%         if RmaskB(i,j)==1 % sample background
%             Bind= Bind+1 ;
%             SampleB(1,Bind)= IColor(i,j);
%         end
% 
%     end
% end
SampleB = IColor(RmaskB)' ; 
Bind = length(SampleB) ; 
%  1-2) compute histogram of F&B pixels--------------------

if Bind >0
    SampleBHist (1,:)= imhist(uint8(SampleB(1,:))) ; % histogram of background samples
    if length(varargin)>0
    SmoothSize =     varargin{1} ;
    else
    SmoothSize=5 ;    
    end
    
    SampleBHistSmooth=SampleBHist ;
    %   1-2-1 ) smoothing the F&B histograms-------------------
    for iter=1 : 1
        for i=1 : 256
            lb= max(1,i-SmoothSize);
            hb= min(256,i+SmoothSize);
            SampleBHistSmooth (1,i)=mean(SampleBHistSmooth(1,lb:hb));
        end
    end
    
    % 1-3) find number of peaks--------------------------------
    
    PeakArea=SmoothSize ;
    SampleBHistPeak=SampleBHistSmooth ;
    PeakNumB=0 ;
    
    for i=1 : 256
        lb= max(1,i-PeakArea);
        hb= min(256,i+PeakArea);
        if (SampleBHistPeak (1,i)== max(SampleBHistSmooth(1,lb:hb)))&&  (SampleBHistPeak (1,i)>0)
            SampleBHistPeak(1,i)= SampleBHistSmooth(1,i) ;
            PeakNumB=PeakNumB+1 ;
        else
            SampleBHistPeak(1,i)= 0 ;
        end
    end
    
end % eof bind==0

PeakNumBSet(channeli) = PeakNumB ;
HistSet = [HistSet;SampleBHist] ; 

end

%% return the result
RPeakNumB=PeakNumBSet ;
RHist = HistSet ; 
% RHist = SampleBHistSmooth ; 


end % eof function