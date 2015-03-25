function [optBoundary poissonMsk] = PoissonOptimalBoundarySparse(imDiff, bwMask)
% SPARSE VERSION OF POISSONOPTIMALBOUNDARYV2.M TO ALLOW FOR LARGER MATRICES AND
% FASTER OPERATION; ALSO STARTS COSTMAT FROM ADJ_8 INSTEAD OF ADJ_4, SINCE
% DESPITE THE LONGER PROCESSING TIME IT PRODUCES LOWER MINCOST FINAL PATH
% IN THE END;
% function to obtain the optimal boundary for poisson editing per the paper
% "Drag and Drop Pasting" by Jia et al;
% Inputs:
%    imDiff: grayscale (fTarget - fSource) where both images are size of cropped
%      source region
%    bwMask: binary mask that segments the object of interest in fSource
% Outputs:
%    optBoundary: the optimal boundary found as list of node numbers
%      (numbering of nodes is columns wise wrt to imDiff matrix, same as matlab linear indexing format)
%    poissonMsk: binary mask where the area within the optimal boundary is
%       set one; to be used to determine the area where poisson equation
%       should be solved; same size as imDiff
% Example: see test_optimal_boundary.m
%DIFFERENCE WITH PREV VERSION IS THAT INPUT IS FIRST NORMALIZED IN RANGE,
%NO NEED TO RESCALE ANY OF THE OUTPUTS BACK IN RANGE SINCE ONE IS BINARY
%AND THE OTHER A RANGE OF INDICES.
tic
graphSearchAlgo = 'Dijkstra_gaimc';

imDiff = double(imDiff);
minvalue = min(imDiff(:));
maxvalue = max(imDiff(:));
imDiff = imDiff - minvalue;
imDiff = imDiff / (maxvalue - minvalue);

connectType = 4; %same as paper, allow 4-connected graph to find optimal boundary
itrMax = 1000;
tol = 0.05; %if percent change in energy is less than tol percent, terminate optimization

% imDiff = fTarget - fSource;

[r,c] = size(imDiff);      
adj_4 = AdjMatrixSparse([r c], connectType);
adj_8 = AdjMatrixSparse([r c], 8);

initialBoundarySub(:,1) =[1:r ones(1,c-2)*r r:-1:1 ones(1, c-2) ]'; %initial boundary is points along the border of the image
initialBoundarySub(:,2) = [ones(1, r) 2:c-1 ones(1,r)*c c-1:-1:2 ]';
initialBoundaryInd = sub2ind([r c], initialBoundarySub(:,1), initialBoundarySub(:,2)); % convert initial boundary to node numbers


currentBoundary = initialBoundaryInd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Preprocessing:check /correct if multiple connected components present in
%%%%% image and try to connect them together; otherwise dijkstra might cut in
%%%%% between the objects thereby eliminating one or more! if multiple
%%%%% objects are found within the mask, it will use the area within convex hull of the objects as the mask 
[L numL] = bwlabel(bwMask);
if numL>1
        L(L>1) = 1; %change all labels to same thing, so that only one convexImage will be returned by regionprops for all objects
        stats = regionprops(L, 'ConvexImage', 'BoundingBox');
        bboxConvHull = round(stats.BoundingBox); %NOTE: coordinates are in matlab image coordinates and should be flipped for matrix indexing!
        convImSize = size(stats.ConvexImage); %stats.ConvexImage is the size of bbox of image
        bwMask(bboxConvHull(2): bboxConvHull(2) + convImSize(1) - 1, bboxConvHull(1): bboxConvHull(1) + convImSize(2) - 1) = stats.ConvexImage; %flipped bbox ul corner since it was in image coordinates (see comment above!)
        disp('--poissonOptimalBoundary: Multiple objects on mask image, using the enclosing convex hull of all objects as mask instead!!!!');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[cutNodes, nodes_cw, nodes_countercw]= FindcutSparse(bwMask); % node lists include the initial boundary point but not the end point located on the mask
bwCutMask = bwMask; bwCutMask(cutNodes) = 1; %cutNodes set to 1 so that costMat is built properly to avoid crossing cutNodes

j = 1;
infCost = 10^5; %arbitrary big number to block nodes that shouldn't be traversed (i.e. nodes on object mask)
K{1} = (1/length(initialBoundaryInd)) * sum(imDiff(initialBoundaryInd)); %1st iteration k initialized according to  initial boundary equal to full boundary
Energy{1} = sum((imDiff(initialBoundaryInd) - K{1}).^2);
energyDiff{1} = tol + 1; % initial value is abitrary number greater than tol
minPath{1} = initialBoundaryInd;
[row_adj col_adj] = find(adj_4);

SID = nodes_cw;
FID = cutNodes;
costs = infCost * ones(length(SID), length(FID)); %initialize to large value
paths = cell(length(SID), length(FID)); %initialize to empty cell array
cutNeighbs = zeros(length(SID), length(FID)); %matrix that shows which SID and FID nodes are 4-neighbors; dijkstra will only run for these combinations of SID/FID, so less number of combos than if 8-neighbors are used
for i = 1:length(SID)
        for z = 1:length(FID)
                if adj_4(SID(i),FID(z))
                        cutNeighbs(i,z) = 1;
                end
        end
end
costsRedux = zeros(size(cutNeighbs)); %initialize

[row_msk col_msk] = find(bwMask);
ind_bwCutMask = find(bwCutMask); %linear indices of nodes on the bwCutMask
% ind_bwCutMask_not = setdiff(1:numel(imDiff), ind_bwCutMask); %linear indices of nodes not on the bwCutMask

if isempty(FID) || length(cutNodes)<2 || min(row_msk(:))==1 || max(row_msk(:)) == size(bwMask,1) || min(col_msk(:))==1 || max(col_msk(:)) == size(bwMask,2)
% if isempty(FID) || length(cutNodes)<2 
        %         error('Selected region is ill-conditioned for optimal boundary calculation; Select regular Poisson paste option');
        disp('-----poissonOptimalBoundary: Selected region is ill-conditioned for optimal boundary calculation; Using regular Poisson paste option!');
        optBoundary = initialBoundaryInd;
        poissonMsk = zeros(r, c);
        poissonMsk(optBoundary) = 1;
        poissonMsk = bwConvexHull(poissonMsk); %obtain mask with pixels within original boundary set to one
%         figure; imshow(poissonMsk);
        
else

        j = 2; %1st iteration handled earlier in the code
        while (j < itrMax) && (abs(energyDiff{j-1}) > tol)
                currentK = K{j - 1};               
                %costMat is defined over three steps in each iteration 
                costMat = adj_4; %initialize; note that in DijkstraSparse a cost of 0 indicates no connections between the nodes; this is unlike dijkstra.m where 0 means cost of 0
                currentDiff = (imDiff(:) - currentK).^2; %need to recompute this in each iteration based on new value of currentK
                costMat = costMat * diag(sparse(currentDiff(:))); %step1: set the ones within each column of costMat (corresponding to j'th node) to the value of currentDiff for that node
                bwcmVector = sparse(1, ind_bwCutMask, ones(length(ind_bwCutMask),1), 1, numel(imDiff)); %sparse vector with a one for each node number that is on bwCutMask
                costMat = costMat * diag(~bwcmVector); %step2: set costs at each column of costMat (corresponding to j'th node) to 0 for all nodes on bwCutMask (i.e. cut their connections), & keep the rest the same (multiply by 1)
                countercwList = []; %initialize
                for i = 1:length(nodes_countercw) %construct matrix of countercw/cutNodes & associated cost of moving if they are 8-neighbors for use in step3 below
                        for z = 1:length(cutNodes)
                                currentSrcNode = nodes_countercw(i);
                                currentTgtNode = cutNodes(z);
                                if adj_4(currentSrcNode, currentTgtNode) == 1
                                        countercwList = [countercwList; currentSrcNode currentTgtNode (imDiff(currentTgtNode)-currentK)^2];
                                end
                        end
                end                
                costMat = costMat + sparse(countercwList(:,1) ,countercwList(:,2) , countercwList(:,3) , size(costMat,1), size(costMat,2)); %step3: cost of going from countercw to cutNodes should be set according to imDiff as well

                for i = 1:length(SID)
                        for z = 1:length(FID)
                                currentSID = SID(i);
                                currentFID = FID(z);
                                if cutNeighbs(i,z) == 1 %compute path only if currentSID & currentFID are 4-neighbors so that path will be a closed loop, and also to save time by reducing computed combinations of SID/FID
                                        if strcmp(graphSearchAlgo, 'DijkstraSparse')
                                                [costs(i,z) paths{i,z}] = DijkstraSparse(costMat,currentSID, currentFID);
                                        elseif strcmp(graphSearchAlgo, 'Dijkstra_gaimc') %this is many times faster than DijkstraSparse, uses a heap in computing shortest path
                                                [currentD currentPred] = Dijkstra_gaimc(costMat, currentSID); %the returned costs are for shortest paths starting from SID to all other nodes in the graph
                                                costs(i,z) = currentD(currentFID); %extract cost of going from currentSID to currentFID;see comment from above line;
                                                paths{i,z} =[]; u = currentFID; while (u ~= currentSID) paths{i,z}=[u paths{i,z}]; u=currentPred(u); end %construct path from list of predecessor nodes returned by the function
                                                paths{i,z} = [u paths{i,z}]; %add the currentSID to the beginning of the path; the condition within while loop above does not let currentSID to be added to path.
                                        end
                                end
                        end
                end                
                sidCosts = (imDiff(SID') - currentK).^2; %costs doesn't include initial cost from SID nodes, so they need to be added
                costsRedux = costs + repmat(sidCosts, [1 length(FID)]);% first add the initial costs
                [minCost minInd] = min(costsRedux(:)); %find min cost/path among all combination of SID/FID points
                minPath{j} = paths{minInd};

                currentBoundary = minPath{j};
                Energy{j} = sum((imDiff(currentBoundary) - currentK).^2);
                energyDiff{j} = 100 * (Energy{j} - Energy{j - 1})/Energy{j-1}; %percent energy change compared to previous best boundary
                K{j} = (1/length(currentBoundary)) * sum(imDiff(currentBoundary));

                j = j + 1;
        end

        optBoundary = minPath{length(minPath)};

        poissonMsk = zeros(r, c);
        poissonMsk(minPath{length(minPath)}) = 1; %actually these will be already included in returned value for inpolygon
%         poissonMskTemp = bwConvexHull(poissonMsk); %Obsolete; obtains mask with pixels within convex hull optimal boundary set to one(convex hull might include undesirable parts!)
        [Xall Yall] = ind2sub([r c], 1:r*c); %find all x,y pairs in the image
        [xv yv] = ind2sub([r c], minPath{length(minPath)}); %the shortest path nodes correspond to the vertices of the polygon that encloses area of interest
        ind_mask = inpolygon(Xall, Yall, xv, yv); %returns logical vector showing whether each of the points in the image are inside/on the polygon or not        
        poissonMsk(ind_mask) = 1;
%         figure(321); subplot(1,3,1); imshow(poissonMskTemp); subplot(1,3,2); imshow(poissonMsk); subplot(1,3,3); asd = zeros(r,c); asd(minPath{length(minPath)})=1;imshow(asd);pause;
end
disp(['Optimal Boundary Found in ' num2str(j-1) ' Iterations; Elapsed Time: ' num2str(toc)]);

%---------------------------------------------------
% Use below to make a video of the optimal boundary through the consecutive
% iterations
%---------------------------------------------------
% h66=figure(66); set(h66, 'Position', [100, 100, 500, 500]);set(h66, 'Color', [1 1 1]); 
% set(h66, 'Color', [1 1 1]);
% imDiffNorm = imDiff - min(imDiff(:));
% imDiffNorm = imDiffNorm / (max(imDiff(:)) - min(imDiff(:)));
% [pathX pathY] = ind2sub([r c],cutNodes);
% [pathcwX pathcwY] = ind2sub([r c],nodes_cw);
% [pathcountercwX pathcountercwY] = ind2sub([r c],nodes_countercw);
% img_temp(:,:,1) = imDiffNorm; img_temp(:,:,2) = imDiffNorm; img_temp(:,:,3) = imDiffNorm; 
% for j=1:length(pathX) %the cut
%         img_temp(pathX(j), pathY(j),1) = 1;
%         img_temp(pathX(j), pathY(j),2) = 0;
%         img_temp(pathX(j), pathY(j),3) = 0;
% end
% for j=1:length(pathcwX) %the one below cut
%         img_temp(pathcwX(j), pathcwY(j),1) = 0;
%         img_temp(pathcwX(j), pathcwY(j),2) = 0;
%         img_temp(pathcwX(j), pathcwY(j),3) = 1;
% end
% imshow(img_temp,[],'InitialMagnification', 'fit'); hold on;
% contor=1;
% mov(contor) = getframe(h66);
% for i = 1:length(minPath)           
%         pause(1);
%         [pathOptX pathOptY] = ind2sub([r c],minPath{i});
% 
%         img_temp(:,:,1) = imDiffNorm; img_temp(:,:,2) = imDiffNorm; img_temp(:,:,3) = imDiffNorm;        
%         for j=1:length(pathOptX)
%                 img_temp(pathOptX(j), pathOptY(j),1) = 1;
%                 img_temp(pathOptX(j), pathOptY(j),2) = 1;
%                 img_temp(pathOptX(j), pathOptY(j),3) = 0;
%         end
%         imshow(img_temp,[],'InitialMagnification', 'fit');
%         title(['Optimal Boundary Iteration = ', num2str(i)],'FontSize',15);    
%         contor = contor+1;
%         mov(contor) = getframe(h66);        
% end
% hold off;
% movie2avi(mov,'/raida/apezeshk/test/temp.avi','compression','None', 'fps',1);
%---------------------------------------------------


%---------------------------------------------------
% Below is obsolete. Resulting video comes out ugly
%---------------------------------------------------
% h66=figure(66);
% set(h66, 'Color', [1 1 1]);
% for i = 1:length(minPath)
%         doMovie = 0;
%         subplot(1,2,1);  imshow(imDiff,[]); hold on;
%         [pathX pathY] = ind2sub([r c],cutNodes);
%         [pathcwX pathcwY] = ind2sub([r c],nodes_cw);
%         [pathcountercwX pathcountercwY] = ind2sub([r c],nodes_countercw);
%         plot(pathY, pathX, 'r-', pathcwY, pathcwX,'y-', pathcountercwY,  pathcountercwX, 'b-'); title('fTarget-fSource with Cut');
%         hold off;
%         subplot(1,2,2); imshow(imDiff, []); hold on;
%         [pathOptX pathOptY] = ind2sub([r c],minPath{i});
%         plot(pathOptY, pathOptX, 'g-'); plot(pathY, pathX, 'r-'); hold off; title(['itr = ', num2str(i)]);
%         pause(1);
%         mov(i) = getframe(h66);
% end
% movie2avi(mov, '/raida/apezeshk/test/optBoundaryMovie.avi','compression','None', 'fps',1);
%---------------------------------------------------