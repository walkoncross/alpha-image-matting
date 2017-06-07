%%
% The Demo of following paper is presented here 
%
% E.Shahrian, D.Rajan, B.Price, S.Cohen, "Improving Image Matting using Comprehensive Sampling Sets" , CVPR 2013 . 
% 
% ------------------------------------------------------------------------
% Running Steps -----------------------------------------------------------
% 1) add source folder to your matlab path
% 2) Load color image with trimap
% 3) call EstimatedAlpha = ComprehensiveSamplingMatting(I, Trimap) ; 
% ------------------------------------------------------------------------




% load image and trimap 
I = imread('Doll.png') ; 
Trimap = imread('Trimap.png'); 

% call ComprehensiveSamplingMatting .
EstimatedAlpha = ComprehensiveSamplingMatting(I, Trimap) ; 


% For further assistance please contact the editor ------------------------
% -------------------------------------------------------------------------
%  Editted by : Ehsan shahrian Varnousfaderani ----------------------------
%  Email : ehsa0004@e.ntu.edu.sg  -----------------------------------------
%  Date : 16 August 2013 --------------------------------------------------
%  Version 1   : 1 October 2014 --------------------------------------------
%  All rights are reserved ! ----------------------------------------------
% -------------------------------------------------------------------------