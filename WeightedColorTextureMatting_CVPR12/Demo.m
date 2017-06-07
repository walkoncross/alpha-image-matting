
% The Demo of following paper is presented here 
%
% E.Shahrian, D.Rajan, Weighted color and texture sample selection for image matting , CVPR 2012, 
% 
% ------------------------------------------------------------------------
% Running Steps -----------------------------------------------------------
% 1) add source folder to your matlab path
% 2) call EstimatedAlpha= WeightedColorTextureMatting(ImagePath, TrimapPath)
% ------------------------------------------------------------------------



% call weighted color and texture matting function 
EstimatedAlpha = WeightedColorTextureMatting('Doll.png', 'Trimap2.png') ; 


% For further assistance please contact the editor ------------------------
% -------------------------------------------------------------------------
%  Editted by : Ehsan shahrian Varnousfaderani ----------------------------
%  Email : ehsa0004@e.ntu.edu.sg  -----------------------------------------
%  Date : 05 March 2013 ---------------------------------------------------
% Version 1 : 1st October 2014 --------------------------------------------
% -------------------------------------------------------------------------