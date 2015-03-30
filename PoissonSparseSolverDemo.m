% Demo file to show operation of lesion/image insertion using poisson image
% editing; See LesionBlenderDemoVideo.avi for how to use this script. Input
% images are expected to be grayscale. 
% Written by: Aria Pezeshk; March 2015
sourceFilename = 'cameraman.tif';
targetFilename = 'cameraman.tif';

%-----------------------------------------
%blend options
normFlag = 0; %flag to determine whether to normalize each source and target image to reside between [0,1] or not
gradMixOption = 1; %0/1, determines whether to do gradient mixing or not
maskOption = 'Optimal'; %'Optimal' for optimal boundary, or 'Standard' to use the whole rectangular user selected ROI;
%-----------------------------------------

source = double(imread(sourceFilename));
if normFlag
        minvalue = min(source(:));
        maxvalue = max(source(:));
        sourceNorm = source - minvalue;
        sourceNorm = sourceNorm / (maxvalue - minvalue);
else
        sourceNorm = source;
end
target = double(imread(targetFilename));
if normFlag
        minvalue = min(target(:));
        maxvalue = max(target(:));
        targetNorm = target - minvalue;
        targetNorm = targetNorm / (maxvalue - minvalue);
else
        targetNorm = target;
end
%------------------------------------------------
% Call up the GUI
%------------------------------------------------
outputParams = simpleBlendGUI(sourceNorm, targetNorm);
rectPos = outputParams.rectPos;
pointPos = outputParams.pointPos;
segMaskFullSize = outputParams.segMaskFullSize;

disp(['Selected rectPos: '  num2str(rectPos)]);
disp(['Selected pointPos: '  num2str(pointPos)]);

%------------------------------------------------
% Find copy paste result
%------------------------------------------------
copyPasteImage = targetNorm;
tlPos = [pointPos(1) - round(rectPos(3)/2), pointPos(2) - round(rectPos(4)/2)]; %top left of insertion area, in matlab image coordinates
copyPasteImage(tlPos(2):(tlPos(2)+rectPos(4)), tlPos(1):(tlPos(1)+rectPos(3))) = imcrop(sourceNorm, rectPos); %copy paste source ROI into target 

%------------------------------------------------
% Test the sparse imageBlender function
%------------------------------------------------
segMask = imcrop(segMaskFullSize, rectPos);

%poisson paste options
pasteAlgoParams.pasteAlgo = 'poisson';
pasteAlgoParams.paramValue = 1; %any value other than -10. which is indicator of error (in master GUI)
pasteAlgoParams.rectPos = rectPos;
pasteAlgoParams.concealOption = 0; %check whether option to conceal is on or off
pasteAlgoParams.postAlphaBlendOption = 0; %check whether post alpha blend option on or off
pasteAlgoParams.targetedGradMixOption = 0; %check whether targeted gradient mixing option on or off
pasteAlgoParams.gradMixOption = 0;
pasteAlgoParams.sliceNodNeighb = [];
pasteAlgoParams.backgroundSliceIndex = 1;
pasteAlgoParams.maxRadiologistMsk{1} = segMask;
% pasteAlgoParams.poissonThresh = poissonThresh; %threshold for convergence of poisson reconstruction
% pasteAlgoParams.omega = 1.5; %NO LONGER NEED THISIN THE SPARSE SOLVER MODE
% pasteAlgoParams.maxNumIter = 3500;
pasteAlgoParams.resizeFactor = 1;
pasteAlgoParams.maskOption = maskOption;
pasteAlgoParams.gradCompOption = 'FullSize';
[imBlend blendOutputs]= ImageBlenderSparse(target, source, pasteAlgoParams, pointPos); %
 
figure; set(gcf, 'Position', get(0,'Screensize')); %force maximize figure
subplot(2,2,1); imshow(sourceNorm,[],  'InitialMagnification', 'fit'); title('Source');
subplot(2,2,2); imshow(targetNorm,[],  'InitialMagnification', 'fit'); title('Target');
subplot(2,2,3); imshow(copyPasteImage,[], 'InitialMagnification', 'fit'); title('Copy and Paste');
subplot(2,2,4); imshow(imBlend, [], 'InitialMagnification', 'fit') ;title('Blending Results Using Sparse Solver');
