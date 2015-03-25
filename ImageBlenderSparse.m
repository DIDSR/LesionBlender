function [imBlend blendOutputs] = ImageBlenderSparse(sliceTargetAll, sliceSourceAll, pasteAlgoParams, pointPos)
% DIFFERENCE WITH IMAGEBLENDER_V2 IS THAT SPARSE SOLVER FOR POISSON
% RECONSTRUCTION IS USED (FOR 2D ONLY, A 3D VERSION HAS TO BE WRITTEN), AS WELL AS SPARSE FAST SOLUTION FOR POISSON
% OPTIMAL BOUNDARY FUNCTION;
% function to perform alpha blending or poisson blending according to input
% parameters
% Inputs:
%     sliceTargetAll: matrix containing the target slice or padded slices (in case of 3d paste)   
%     sliceTargetAll: matrix containing the source slice or padded slices (in case of 3d paste)  
%    pasteAlgo: 'alpha' or 'poisson' or '3dpoisson'
%    pointPos: center point of where cropped img should be pasted onto
%          target; assumes image coordinates
%    pasteAlgoParams: structure containing the necessary params for the
%       various paste modes
% Output:
%    imBlend: the result of blending algorithm. in order to adhere to the
%         format in the GUI that is used for window/level sliders to work correctly,
%         image is type cast to int16  
normalizeFlag = 0; %flag to determine whether to normalize each of source and target to be between [0,1] before further processing
borderGradientFlag = 1;

sliceSourceAll = double(sliceSourceAll);
sliceTargetAll = double(sliceTargetAll);

if normalizeFlag==1
        minvalueSrc = min(sliceSourceAll(:)); %normalize source to be in [0 1]
        maxvalueSrc = max(sliceSourceAll(:));
        sliceSourceNorm = sliceSourceAll - minvalueSrc;
        sliceSourceNorm = sliceSourceNorm / (maxvalueSrc - minvalueSrc);

        minvalue = min(sliceTargetAll(:)); % normalize target to be in [0 1]
        maxvalue = max(sliceTargetAll(:));
        SliceTargetNorm = sliceTargetAll - minvalue;
        SliceTargetNorm = SliceTargetNorm / (maxvalue - minvalue);
else
        minvalueSrc = 0; maxvalueSrc = 1;
        minvalue = 0; maxvalue = 1;
        sliceSourceNorm = sliceSourceAll;
        SliceTargetNorm = sliceTargetAll;
end

rectPos = round(pasteAlgoParams.rectPos); % round in case of non-integer coordinates
srcWidth = rectPos(3) + 1; %3rd element + 1 gives width
srcHeight = rectPos(4) + 1; %4th element + 1 gives height

centerX = round(pointPos(2)); %matlab image and matrix coordinates are switched
centerY = round(pointPos(1));

srcTop = rectPos(2); %limits of the source crop box
srcBottom = rectPos(2) + rectPos(4);
srcLeft = rectPos(1);
srcRight = rectPos(1) + rectPos(3);

tgtTop = centerX - floor(srcHeight/2); % coordinates of the perimeter of bounding box of where img will be pasted, according to matrix coordinates (rather than matlab image coordinates)
tgtBottom = tgtTop + srcHeight - 1;
tgtLeft = centerY - floor(srcWidth/2);
tgtRight = tgtLeft + srcWidth - 1; 

tileImage = []; %initialize; if option to conceal is selected, this will contain seamless tile of background; 

pasteAlgo = pasteAlgoParams.pasteAlgo;
paramValue = pasteAlgoParams.paramValue; %for case of 'alpha' should be between [0,1], ; for erroneous input parameters, it will be set to -10 in ctm
maskOption = pasteAlgoParams.maskOption; % 'Standard' or 'Optimal' in case of pasteAlgo==poisson
concealOption = pasteAlgoParams.concealOption; % 0/1 depending on whether conceal option is desired
postAlphaBlendOption = pasteAlgoParams.postAlphaBlendOption;% 0/1 depending on whether postAlphaBlend option is desired
gradMixOption = pasteAlgoParams.gradMixOption; %check whether regular grad mixing option is on or off
targetedGradMixOption = pasteAlgoParams.targetedGradMixOption;%check whether targeted gradient mixing option on or off
backgroundSliceIndex = pasteAlgoParams.backgroundSliceIndex; %index of which src slice to crop the background region from
sliceNodNeighb = pasteAlgoParams.sliceNodNeighb; %vector with same length as number of slices to paste; 1 indicates slice index is adjacent to start/end of nodule
maxRadiologistMsk = pasteAlgoParams.maxRadiologistMsk; %max radiologist mask is same size as cropped area

if isfield(pasteAlgoParams, 'contrastScaleFactor') %if passed to this function, this field defines a scaling of the grad field of source by value of 'contrastScaleFactor', s.t. blending changes contrast of insertion; affects only poisson/2dmultislice
        contrastScaleFactor = pasteAlgoParams.contrastScaleFactor;
else
        contrastScaleFactor = 1;
end
if isfield(pasteAlgoParams, 'resizeFactor') %if passed to this function, this field defines a scaling of the grad field of source by value of 'contrastScaleFactor', s.t. blending changes contrast of insertion; affects only poisson/2dmultislice
        resizeFactor = pasteAlgoParams.resizeFactor; %if set to 0<value<1, image will be first resized by that scale to compute optimal boundary faster; reconstruction done at full original image size;
        if (resizeFactor > 1) || (resizeFactor <= 0)
                error('function imageBlenderSparse: resizeFactor provided in pasteAlgoParams structure should be between 0 and 1');
        end
else
        resizeFactor = 1;
end

boundaryPixDefault(:,1) =[1:srcHeight ones(1,srcWidth-2)*srcHeight srcHeight:-1:1 ones(1, srcWidth-2) ] + tgtTop - 1; %initial boundary is points along the border of the cropped image, so convert that to target coordinates; will be changed if Optimal boundary option selected
boundaryPixDefault(:,2) = [ones(1, srcHeight) 2:srcWidth-1 ones(1,srcHeight)*srcWidth srcWidth-1:-1:2 ] + tgtLeft - 1;
switch pasteAlgo
        case {'Alpha','alpha'}
                if paramValue<0 || paramValue>1
                        disp('-----Error!');
                        imBlend = zeros(size(SliceTargetNorm));
                else
%                         imBlend = double(SliceTargetNorm);
%                         sourceCroppedNorm = imcrop(sliceSourceNorm, rectPos);
%                         imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,:) = paramValue*imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,:) + (1 - paramValue) * sourceCroppedNorm; % note that imBlend is type cast to int16 before being returned by the function

                        imBlend = double(sliceTargetAll); %the commented part does the alpha blending based on non-normalized images; it has to normalize the blend result as in the end of the file there is a non-normalizing factor applied to go back to HU
                        sourceCroppedNorm = imcrop(sliceSourceAll, rectPos);
                        imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,:) = paramValue*imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,:) + (1 - paramValue) * sourceCroppedNorm; % note that imBlend is type cast to int16 before being returned by the function
                        temp = imBlend - minvalue;
                        temp = temp/(maxvalue-minvalue); imBlend=temp;                     
                end
                
                boundaryPixAll = {boundaryPixDefault};

        case {'Poisson','poisson','2dmultislicepoisson','2DMultiSlicePoisson'}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%% Insertion into SliceTarget Image Using Poisson Blend
                %%%%%%%%%% Handles both single/multiSlice cases
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if paramValue==-10 %in ctm, when invalid src/tgt/background slices are selected in multislice 2d poisson, paramValue is set to -10
                        disp('-----Error!');
                        imBlend = zeros(size(SliceTargetNorm));
                        boundaryPixAll = cell(1, size(SliceTargetNorm,3)); %initialize to the default boundaryPix
                        boundaryPixAll(:) = {boundaryPixDefault};
                else
                        numSlicesToPaste = size(SliceTargetNorm, 3);
                        boundaryPixAll = cell(1, numSlicesToPaste);
                        imBlend = SliceTargetNorm; %preallocate imBlend
                        sourceCroppedNorm = sliceSourceNorm(srcTop: srcBottom, srcLeft:srcRight, :); %if concealOption==1, this will be modified below

                        %%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%Preprocessing
                        if concealOption %if conceal option is activated, conceal attachments/artifacts in cropped region; both sliceSourceNorm/sourceCroppedNorm will be modified with updated values!!!
%                                 maxRadiologistMsk = pasteAlgoParams.maxRadiologistMsk; %max radiologist mask is same size as cropped area
                                rectPosBackground  = round(pasteAlgoParams.rectPosBackground); % round in case of non-integer coordinates
                                bkgrndTop = rectPosBackground(2);
                                bkgrndBottom = rectPosBackground(2) + rectPosBackground(4);
                                bkgrndLeft = rectPosBackground(1);
                                bkgrndRight = rectPosBackground(1) + rectPosBackground(3);
                                backgroundCroppedNorm = sliceSourceNorm(bkgrndTop:bkgrndBottom, bkgrndLeft:bkgrndRight, backgroundSliceIndex); %crop the bgnd region from correct slice
                                concealOpMode = 'poissonReplacePSF'; %'poissonReplacePointwisePSF', 'poissonReplacePSF','poissonReplace', 'directReplace', 'membraneInterpolant', 'exemplarInpainting'
                                medFiltFlag = 1;
                                [sourceCroppedNorm tileImage]= concealerV2(sourceCroppedNorm, backgroundCroppedNorm, maxRadiologistMsk, concealOpMode, medFiltFlag);
                                sliceSourceNorm(srcTop: srcBottom, srcLeft:srcRight, :) = sourceCroppedNorm; %NOTE: sliceSourceNorm being updated here with the cropped concealed content
                                if borderGradientFlag==1 %if flag is set, set the following row/column to mirror the inside of ROI, s.t. gradient will be zero to avoid post conceal issues across the boundary
                                        sliceSourceNorm(srcBottom+1, srcLeft:srcRight, :) = sliceSourceNorm(srcBottom, srcLeft:srcRight, :);
                                        sliceSourceNorm(srcTop: srcBottom, srcRight+1, :) = sliceSourceNorm(srcTop: srcBottom, srcRight, :); 
                                end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%

                        for s = 1:numSlicesToPaste %do slice by slice 2D poisson blending

                                boundaryPix = boundaryPixDefault;
                                % currentSourceCroppedNorm = imcrop(sliceSourceNorm(:,:,s), rectPos);
                                currentSourceCroppedNorm = sourceCroppedNorm(:,:,s);
                                [TgtH TgtV] = imgrad(SliceTargetNorm(:,:,s)); % find the img gradients
                                [SrcHFull SrcVFull] = imgrad(contrastScaleFactor*sliceSourceNorm(:,:,s)); %to avoid boundary issues, compute gradient from full source then crop; sliceSourceNorm is updated above with the cropped concealed content
                                %ACHTUNG:multiplication of grad field by
                                %contrastScaleFactor above changes the contrast
                                %of insertion during poisson recon. 
                                SrcH = imcrop(SrcHFull, rectPos); SrcV = imcrop(SrcVFull, rectPos);
%                                 [SrcH SrcV] = imgrad(currentSourceCroppedNorm);

                                currentImBlend = SliceTargetNorm(:,:,s); % initialize currentImBlend, it will later become the copy with sourceCropped region directly pasted into target
                                Th = TgtH;
                                Tv = TgtV;

                                if strcmpi(maskOption, 'Standard') %case of standard poisson where whole area is reconstructed
                                        msk = zeros(size(currentImBlend));
                                        msk(tgtTop:tgtBottom,tgtLeft:tgtRight,:) = 1;
                                        optMsk = ones(size(currentSourceCroppedNorm)); %

                                elseif strcmpi(maskOption, 'Optimal') %
                                        fSource = currentSourceCroppedNorm;
                                        fTarget = SliceTargetNorm(tgtTop:tgtBottom, tgtLeft:tgtRight, s);
                                        imDiff = fTarget - fSource;

                                        if isempty(pasteAlgoParams.maxRadiologistMsk{s}) %if maximal radiologist marking for current slice available, use that; otherwise use original bbox  to determine cost matrix weighting(can replace with rough segmentation)
                                                %level =graythresh(fSource);
                                                % bwMask = im2bw(fSource, level); %rough segmentation of object of interest, prob should be replaced with sth like GrabCut; to be used in building the adjacency matrix
                                                msk = zeros(size(currentImBlend));
                                                msk(tgtTop:tgtBottom,tgtLeft:tgtRight,:) = 1;
                                                optMsk = ones(size(currentSourceCroppedNorm)); %
                                        else
                                                bwMask = pasteAlgoParams.maxRadiologistMsk{s};
                                                [r_imDiff c_imDiff] = size(imDiff); %imDiff & bwMask are same size as original selected ROI

                                                if resizeFactor==1 %if operating on original image size for optimal boundary, then just use original image
                                                        [optNodes optMsk] = PoissonOptimalBoundarySparse(imDiff, bwMask);
                                                        optNodesSub = bwboundaries(optMsk,4); %can also use optNodes directly, but for conformity will do same thing as when we have to resize the image below
                                                        optMsk = double(optMsk);
                                                        if length(optNodesSub)>1 %variable is a cell array for indices of boundaries of all objects; in case multple objects are found, that is an error
                                                                error('function ImageBlenderSparse: multiple objects being found within optimal boundary, something is wrong!');
                                                        end
                                                        optNodes = sub2ind([r_imDiff c_imDiff], optNodesSub{1}(:,1), optNodesSub{1}(:,2)); %1st column is row numbers, 2nd column is col numbers; convert to linear indices.
                                                else %if operating on smaller image size for optimal boundary, additional steps should be taken                                                        
                                                        imDiffRsz = imresize(imDiff, resizeFactor);
                                                        bwMaskRsz = double(imresize(logical(bwMask), resizeFactor)); %logical is done so that values of resized image are still 0 and 1
                                                        [optNodesRsz optMskRsz] = PoissonOptimalBoundarySparse(imDiffRsz, bwMaskRsz); %compute optimal boundary off of the resized inputs
                                                        optMsk = imresize(logical(optMskRsz), [r_imDiff c_imDiff]); %revert optMsk to same size as imDiff
                                                        optNodesSub = bwboundaries(optMsk,4); %optNodesRsz cannot be converted to get the list of nodes, so find them from optMsk using bwboundaries
                                                        optMsk = double(optMsk);
                                                        if length(optNodesSub)>1 %variable is a cell array for indices of boundaries of all objects; in case multple objects are found, that is an error
                                                                error('function ImageBlenderSparse: multiple objects being found within optimal boundary, something is wrong!');
                                                        end
                                                        optNodes = sub2ind([r_imDiff c_imDiff], optNodesSub{1}(:,1), optNodesSub{1}(:,2)); %1st column is row numbers, 2nd column is col numbers; convert to linear indices.
                                                end
                                                % optimalBoundary = minPath{length(minPath)};
                                                [optNodesSubX optNodesSubY] = ind2sub([srcHeight srcWidth], optNodes); %node numbers returned as boundary are wrt cropped image size, so next convert them into actual full target coordinates
                                                optNodesSubX = optNodesSubX + tgtTop - 1;
                                                optNodesSubY = optNodesSubY + tgtLeft - 1;
                                                boundaryPix = []; boundaryPix(:,1) = optNodesSubX; boundaryPix(:,2) = optNodesSubY;
                                                
                                                msk = zeros(size(currentImBlend));
                                                msk(tgtTop:tgtBottom,tgtLeft:tgtRight,:) = optMsk; %msk is same size as whole target image, optMsk is same size as cropped source
                                        end
                                end

                                msk_logical = logical(msk); optMsk_logical = logical(optMsk); %logical versions of variables so that they can be used for matrix indexing
                                currentImBlend(msk_logical) = currentSourceCroppedNorm(optMsk_logical);  %only replace target with parts of source that are within boundary
                               
                                % in regular mode, match gradient of area within target area  to that of source
                                Th(msk_logical)  = SrcH(optMsk_logical);  % initialize; in regular mode, match gradient of area within target area  to that of source
                                Tv(msk_logical) = SrcV(optMsk_logical); % initialize; in regular mode, match gradient of area within target area  to that of source
                                        
                                if gradMixOption || targetedGradMixOption % for gradient mixing, match gradient of area within target to target/source dependingon which has larger abs value
                                        Th(msk_logical)  = SrcH(optMsk_logical); %first replace gradient of target with source within designated area
                                        Tv(msk_logical) = SrcV(optMsk_logical);
                                        srcOrTgt = (TgtH.^2 + TgtV.^2)> (Th.^2 + Tv.^2); %now compare gradient values of original target vs. new target for all points (we aren't changing gradient values outside of designated area, so comparison against full image is ok)                                        
                                        if targetedGradMixOption && ~isempty(maxRadiologistMsk{s}) && ~sliceNodNeighb(s)%in targeted gradient mixing areas on radiologist mask maintain src gradient (so that veins from tgt passing thru nodule don't get artificially super bright)
                                                srcOrTgt(tgtTop:tgtBottom,tgtLeft:tgtRight) = srcOrTgt(tgtTop:tgtBottom,tgtLeft:tgtRight) & ~maxRadiologistMsk{s};
                                        end                                        
                                        srcOrTgt = logical(srcOrTgt);
                                        Th(srcOrTgt) = TgtH(srcOrTgt); %now replace gradient within designated area with max of gradient of target/source
                                        Tv(srcOrTgt) = TgtV(srcOrTgt);
                                end

                                Y2 = PoissonSparseSolver(currentImBlend, Th, Tv, msk);

                                if postAlphaBlendOption && concealOption %if postBlend is selected, blend more target background into source background
                                        lowerBlendValue = 0.7; upperBlendValue = 0.9;
                                        currentSliceTargetNormCropped = SliceTargetNorm(tgtTop:tgtBottom,tgtLeft:tgtRight,s);
                                        Y2Cropped = Y2(tgtTop:tgtBottom,tgtLeft:tgtRight);
                                        currentMaxRadMsk = maxRadiologistMsk{s};
                                        temp = Y2(tgtTop:tgtBottom,tgtLeft:tgtRight); %initialize
                                        imBlend(:,:,s) = Y2; % note that imBlend is type cast to int16 before being returned by the function                                        

                                        temp(~currentMaxRadMsk) = lowerBlendValue * Y2Cropped(~currentMaxRadMsk) + (1-lowerBlendValue) * currentSliceTargetNormCropped(~currentMaxRadMsk);
                                        temp(currentMaxRadMsk) = upperBlendValue * Y2Cropped(currentMaxRadMsk) + (1-upperBlendValue) * currentSliceTargetNormCropped(currentMaxRadMsk);

                                        imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,s) = temp;
                                else %otherwise uniformly do a minor alpha blend over entire paste region
                                        upperBlendValue = 1;
                                        imBlend(:,:,s) = Y2; % note that imBlend is type cast to int16 before being returned by the function                                        
                                        imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,s) = upperBlendValue * Y2(tgtTop:tgtBottom,tgtLeft:tgtRight) + (1 - upperBlendValue) * SliceTargetNorm(tgtTop:tgtBottom,tgtLeft:tgtRight,s); %optionally end result could be alpha blended for additional smoothing
                                end
                                
                                boundaryPixAll{s} = boundaryPix;                                
                        end
                end

              
                
        case {'3DPoisson','3dpoisson'} %in case of 3d paste, SliceTargetNorm/SliceSourceNorm will have several pages each (they are 3d matrices)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%% Insertion into SliceTarget Image Using 3DPoisson Blend
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%NOTE: USING sliceNodNeighb  TO GUIDE HOW TARGETED GRAD
                %%MIXING IS PERFORMED HAS NOT BEEN IMPLEMENTED IN 3D
                %%YET!!!!!!! ALSO A 5% ALPHA BLENDING IS CURRENTLY BEING
                %%DONE AT THE END OF 3DPOISSON BY DEFAULT!!!!!
                
                if paramValue==-10 %in ctm, when invalid src/tgt/background slices are selected paramValue is set to -10 
                        disp('-----Error!');
                        imBlend = zeros(size(SliceTargetNorm));
                        boundaryPixAll = cell(1, size(SliceTargetNorm,3)); %initialize to the default boundaryPix
                        boundaryPixAll(:) = {boundaryPixDefault};
                else
                        imBlend = SliceTargetNorm; % initialize imBlend, it will later become the copy with source region directly pasted into target
                        sourceCroppedNorm = sliceSourceNorm(srcTop: srcBottom, srcLeft:srcRight, :);
                        
                        boundaryPixAll = cell(1, size(SliceTargetNorm,3)); %initialize to the default boundaryPix
                        boundaryPixAll(:) = {boundaryPixDefault};
                        rangeZ = 3:size(SliceTargetNorm,3)-2; %source/target padded with 2 slices on top and bottom each; this will be used to have non-zero masks only on pages that we want to solve for (sans padding)
                        concealRangeZ = rangeZ; %initialize concealRangeZ same as rangeZ; if conceal is on will allow one padding on top/bottom to be included in 3d reconstruction (altho they won't be returned in imBlend))
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%Preprocessing
                        if concealOption %if conceal option is activated, conceal attachments/artifacts in cropped region; both sliceSourceNorm/sourceCroppedNorm will be modified with updated values!!!
                                maxRadiologistMsk = pasteAlgoParams.maxRadiologistMsk;%max radiologist mask is same size as cropped area
                                rectPosBackground  = round(pasteAlgoParams.rectPosBackground); % round in case of non-integer coordinates
                                bkgrndTop = rectPosBackground(2);
                                bkgrndBottom = rectPosBackground(2) + rectPosBackground(4);
                                bkgrndLeft = rectPosBackground(1);
                                bkgrndRight = rectPosBackground(1) + rectPosBackground(3);
                                backgroundCroppedNorm = sliceSourceNorm(bkgrndTop:bkgrndBottom, bkgrndLeft:bkgrndRight, backgroundSliceIndex);
                                concealOpMode = 'poissonReplacePSF'; %'poissonReplacePSF', 'poissonReplace', 'directReplace', 'membraneInterpolant'
                                medFiltFlag = 1;
                                [sourceCroppedNorm tileImage]= concealer3D(sourceCroppedNorm, backgroundCroppedNorm, maxRadiologistMsk, rangeZ, concealOpMode, medFiltFlag);
                                sliceSourceNorm(srcTop: srcBottom, srcLeft:srcRight, :) = sourceCroppedNorm;
                                concealRangeZ = (rangeZ(1)-1): (rangeZ(end)+1); %adds one more slice on each stack to be reconstructed for improved boundary conditions in z direction
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%

                        physicalSpacing = [1 1 1];  %physicalSpacing = physicalSpacing/sum(physicalSpacing);
                        [TgtX TgtY TgtZ] = imgrad3D(SliceTargetNorm, physicalSpacing);
                        Tx = TgtX; Ty = TgtY; Tz = TgtZ;  
                        auxSlices = setdiff(concealRangeZ, rangeZ); %index of "smooth encasing" slices that should be included in reconstruction (but won't be returned in imBlend)
                        
                        if ~concealOption %in regular mode, crop gradients from grad of full size image stack
                                [SrcXFull SrcYFull SrcZFull] = imgrad3D(sliceSourceNorm, physicalSpacing);
                                SrcX = SrcXFull(srcTop: srcBottom, srcLeft:srcRight, :);
                                SrcY = SrcYFull(srcTop: srcBottom, srcLeft:srcRight, :);
                                SrcZ = SrcZFull(srcTop: srcBottom, srcLeft:srcRight, :);
                        else %in conceal mode boundaries might have altered wrt full size, so ompute src grad from the cropped image stack
                                [SrcX SrcY SrcZ] = imgrad3D(sourceCroppedNorm, physicalSpacing);
                        end

                        %                         msk = zeros(size(imBlend));
                        %                         msk(tgtTop:tgtBottom,tgtLeft:tgtRight,rangeZ) = 1;
                        %                         msk_logical = logical(msk);
                        %                         optMsk_logical = true(size(sourceCroppedNorm));
                        %                         imBlend(msk_logical) = sourceCroppedNorm(optMsk_logical);  %only replace target with parts of source that are within boundary

                        if strcmpi(maskOption, 'Standard') %case of standard poisson where whole area in each plane is reconstructed
                                msk = zeros(size(imBlend));
                                msk(tgtTop:tgtBottom,tgtLeft:tgtRight,concealRangeZ) = 1;
                                optMsk = zeros(size(sourceCroppedNorm));
                                optMsk(:,:, concealRangeZ) = repmat(ones(srcHeight, srcWidth), [1 1 numel(concealRangeZ)]);
                               
                        elseif strcmpi(maskOption, 'Optimal') %
                                optMsk = zeros(size(sourceCroppedNorm)); %initialize
                                msk = zeros(size(imBlend)); %initialize
                                
                                if concealOption %in case of conceal, 1st top/bottom paddin slices should also be reconstructed, but won't need to go thru optimal boundary computation                                        
                                        optMsk(:,:,auxSlices) = repmat(ones(srcHeight, srcWidth), [1 1 numel(auxSlices)]);
                                        msk(tgtTop:tgtBottom,tgtLeft:tgtRight,auxSlices) = 1;
                                end
                                
                                for s = rangeZ %loop thru primary slices that need to be solved for, and for each find the optMsk and optimal Boundary
                                        boundaryPix = boundaryPixDefault; %initialize
                                        fSource = sourceCroppedNorm(:,:,s);
                                        fTarget = SliceTargetNorm(tgtTop:tgtBottom, tgtLeft:tgtRight, s);
                                        imDiff = fTarget - fSource;

                                        if isempty(pasteAlgoParams.maxRadiologistMsk{s}) %if maximal radiologist marking for current slice available, use that; otherwise use original bbox  to determine cost matrix weighting(can replace with rough segmentation)
                                                %level =graythresh(fSource);
                                                % bwMask = im2bw(fSource, level); %rough segmentation of object of interest, prob should be replaced with sth like GrabCut; to be used in building the adjacency matrix
                                                msk(:,:,s) = zeros(size(imBlend(:,:,s)));
                                                msk(tgtTop:tgtBottom,tgtLeft:tgtRight,s) = 1;
                                                optMsk(:,:,s) = ones(size(sourceCroppedNorm(:,:,s))); %
                                        else
                                                bwMask = pasteAlgoParams.maxRadiologistMsk{s};

                                                [optNodes currentOptMsk] = PoissonOptimalBoundarySparse(imDiff, bwMask);
                                                % optimalBoundary = minPath{length(minPath)};
                                                [optNodesSubX optNodesSubY] = ind2sub([srcHeight srcWidth], optNodes); %node numbers returned as boundary are wrt cropped image size, so next convert them into actual full target coordinates
                                                optNodesSubX = optNodesSubX + tgtTop - 1;
                                                optNodesSubY = optNodesSubY + tgtLeft - 1;
                                                boundaryPix = []; boundaryPix(:,1) = optNodesSubX; boundaryPix(:,2) = optNodesSubY;

                                                msk(:,:,s) = zeros(size(imBlend(:,:,s)));
                                                msk(tgtTop:tgtBottom,tgtLeft:tgtRight,s) = currentOptMsk; %msk is same size as whole target image, optMsk is same size as cropped source
                                                optMsk(:,:,s) = currentOptMsk;
                                        end
                                        boundaryPixAll{s} = boundaryPix;
                                end
                        end
                        msk_logical = logical(msk); optMsk_logical = logical(optMsk); %logical versions of variables so that they can be used for matrix indexing
                        imBlend(msk_logical) = sourceCroppedNorm(optMsk_logical);  %only replace target with parts of source that are within boundary
                                      

                        % in regular mode, match gradient of area within target area  to that of source
                        Tx(msk_logical)  = SrcX(optMsk_logical); % initialize; in regular mode, match gradient of area within target area  to that of source
                        Ty(msk_logical) = SrcY(optMsk_logical);
                        Tz(msk_logical) = SrcZ(optMsk_logical);
                        if gradMixOption || targetedGradMixOption  % for gradient mixing, match gradient of area within target to target/source dependingon which has larger abs value
                                Tx(msk_logical)  = SrcX(optMsk_logical); %first replace gradient of target with source within designated area
                                Ty(msk_logical) = SrcY(optMsk_logical);
                                Tz(msk_logical) = SrcZ(optMsk_logical);
                                srcOrTgt = (TgtX.^2 + TgtY.^2 + TgtZ.^2)> (Tx.^2 + Ty.^2 + Tz.^2); %now compare gradient values of original target vs. new target for all points (we aren't changing gradient values outside of designated area, so comparison against full image is ok)
                                if targetedGradMixOption %in targeted gradient mixing areas on radiologist mask maintain src gradient (so that veins from tgt passing thru nodule don't get artificially super bright)
                                        for k = 1:size(sourceCroppedNorm,3)     %go thru every slice and compare against the maximal radiologist mask
                                                if ~isempty(maxRadiologistMsk{k}) %padding slices have empty msk
                                                        srcOrTgt(tgtTop:tgtBottom,tgtLeft:tgtRight,k) = srcOrTgt(tgtTop:tgtBottom,tgtLeft:tgtRight,k) & ~maxRadiologistMsk{k};
                                                end
                                        end
                                end
                                srcOrTgt = logical(srcOrTgt);
                                Tx(srcOrTgt) = TgtX(srcOrTgt); %now replace gradient within designated area with max of gradient of target/source
                                Ty(srcOrTgt) = TgtY(srcOrTgt);
                                Tz(srcOrTgt) = TgtZ(srcOrTgt);
                        end
                        
                        Y2 = PoissonGaussSeidel3D( imBlend, Tx, Ty, Tz, msk ,physicalSpacing, 1024, 1e-3, 1.5, true); %in case of 'conceal' option, Y2 will contain reconstructed padding slices, but imBlend will contain only central slices!                        
                        imBlend(:,:,concealRangeZ) = Y2(:,:,concealRangeZ); % note that imBlend is type cast to int16 before being returned by the function; note that original pad slices are returned by this (not the concealed/reconstructed pad slices)
                        if concealOption %in this case the one extra padding slices from conceal reconstruction should stay same as original slices in the target
                                imBlend(:,:, auxSlices) = SliceTargetNorm(:,:,auxSlices);
                        end
                        imBlend(tgtTop:tgtBottom,tgtLeft:tgtRight,rangeZ) = 0.95 * Y2(tgtTop:tgtBottom,tgtLeft:tgtRight, rangeZ) + (1 - 0.95) * SliceTargetNorm(tgtTop:tgtBottom,tgtLeft:tgtRight,rangeZ); %optionally end result could be alpha blended for additional smoothing 
                end
end

% if normalizeFlag==1
%         imBlend = int16(imBlend*(maxvalue - minvalue) + minvalue);
% else
%         imBlend = int16(imBlend);
% end
blendOutputs.boundaryPixAll = boundaryPixAll;
blendOutputs.tileImage = tileImage;
% blendOutputs.tileImage =  int16(tileImage*(maxvalue - minvalue) + minvalue);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Deletion from Source Image & Replacement with Somewhere Else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DX = 359; % coordinates of top left of region we want to paste from source into the part being deleted from it 
% DY = 309;
% 
% holeFill = sourceNorm;
% holeFill(GY:GY+h,GX:GX+w,:) = sourceNorm(DY:DY+h,DX:DX+w,:);
% Th = SrcH; % re-initialize the forward gradient fields to source image
% Tv = TgtV;
% Th(GY:GY+h,GX:GX+w,:) = SrcH(DY:DY+h,DX:DX+w,:); % copying gradient field from different region within same image onto different part of the image
% Tv(GY:GY+h,GX:GX+w,:) = TgtV(DY:DY+h,DX:DX+w,:);
% msk_del = zeros(size(holeFill));
% msk_del(GY:GY+h,GX:GX+w,:) = 1;
% holeFillFinal = PoissonJacobi( holeFill, Th, Tv, msk_del );
% figure(6); subplot(1,2,1); imshow(holeFill,[]); subplot(1,2,2); imshow(holeFillFinal,[]);


