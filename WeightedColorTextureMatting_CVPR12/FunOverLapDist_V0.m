function [ROverLap] = FunOverLapDist_V0 (FSample , BSample, NumBin, BinMinSample_Percent)

% FSamples and BSample should be row wise , this means that samples are in
% the rows and channels are columns.  therefore FSample [20x3] have 20
% samples for three channels .
%
% NumBin = indicate number of bins
%BinMinSample_Percent = indicate minimum number of samples for bins based
%of percent .
%

[numFSample, NumChannel] = size(FSample) ;
if NumChannel >1

    for k=1 : NumChannel
        
        FSample_Channel = FSample(:,k) ; 
        BSample_Channel = BSample(:,k) ; 
        TOverLap(k)= FunOverLapDist_V0 (FSample_Channel , BSample_Channel, NumBin, BinMinSample_Percent) ; 
        
      
        
    end
    % Return The Result 
    ROverLap = mean(TOverLap) ; 
    
else
    
    
    
    FHist = imhist (uint8(FSample), NumBin) ;
    BHist = imhist (uint8(BSample), NumBin) ;
    
    BinMinSamples4F = BinMinSample_Percent *sum(FHist) ;
    BinMinSamples4B = BinMinSample_Percent *sum(BHist) ;
    
    FHist (FHist<BinMinSamples4F)=0 ;
    BHist (BHist<BinMinSamples4B)=0 ;
    
    FHist = FHist / sum(FHist) ;
    BHist = BHist / sum(BHist) ;
    
    
    % OverLap Computation -----------------------------------------------------
    HistInterMask= ((FHist>0 )& (BHist>0 )) ;
    
    Overlapping = (sum(FHist(HistInterMask).*BHist(HistInterMask)))/ (sum(FHist.^2+BHist.^2)/2) ;
    
    
    % Return The Result -------------------------------------------------------
    ROverLap = Overlapping ;
    
end  % EOF IF NumChannels










