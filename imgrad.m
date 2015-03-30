function [forwardH forwardV backwardH backwardV] = imgrad(inputImg)
% function to calculate horizontal and vertical gradients in image.
% Outputs:
%      forwardH: forward horizontal difference
%      forwardV: forward vertical difference
%      backwardH: forward horizontal difference
%      backwardV: forward vertical difference
%
% Inputs:
%      inputImg: input image (grayscale)
kernelH = [ 0,-1, 1 ];
kernelV = [ 0;-1; 1 ];

forwardH = imfilter(inputImg,kernelH,'replicate');
forwardV = imfilter(inputImg,kernelV,'replicate');

if( nargout >= 3 )
 backwardH = circshift(forwardH,[0,1]);
 backwardV = circshift(forwardV,[1,0]);
end
