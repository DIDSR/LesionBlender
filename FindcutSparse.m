function  [cutNodes, nodes_cw, nodes_countercw]= FindcutSparse(bwMask)
% NEW FUNCTION BASED ON FINDCUT.M, WITH THE DIFFERENCE BEING THAT
% IMPLEMENTATION IS NOW BASED ON SPARSE MATRICES FOR SPEEDUP;
% for now use a simple implementation: cut is taken from right hand center
% of image, to nearest point on the mask; function returns the list of
% nodes located on the cut; cutNodes and other outputs initially include the initial and last points
% as well, but last point which is on the mask is dropped in the returned
% values from function
%
% Inputs:
%   bwMask: binary mask where pixels on object of interest are set to one;
%     for now this comes from binraization of image using otsu, and
%     therefore it could contain background pixels, and/or multiple objects in
%     the ROI
% Outputs:
%   cutNodes: nodes located on the cut; for now these are nodes starting
%      from right hand center of image, and ending on closest point on mask
%      (excluding the end node)
%   nodes_cw: nodes just to the side of the cut (clockwise); excludes the
%      end node
%   nodes_countercw: nodes just to the side of the cut (counter-clockwise); excludes the
%      end node
%
% Example:
%   bwMask = zeros(25);
%   bwMask(2:5, 5:18) = 1;
%   nodes  = findcut(bwMask)
% Written By: Aria Pezeshk; 2014

% Later, initial position could be looped around
% the boundary and the cut with shortest overall distance would be picked

[h w] = size(bwMask);
redgeCentX = floor(h/2);
redgeCentY = w;

[row col] = find(bwMask);
redgeCentMat = repmat([redgeCentX redgeCentY], length(row), 1);

distMat = ((row - redgeCentMat(:,1)).^2) + ((col - redgeCentMat(:,2)).^2);
[temp indMinDist] = min(distMat);

adj_4 = AdjMatrixSparse([h w], 4); %get adjacency matrix for 4 or 8-neighbor connectivity 
adj_8 = AdjMatrixSparse([h w], 8);

costMat = adj_4 + sqrt(2) * (adj_8 - adj_4); %move across 4 neighbors with cost of 1, and along 8 neighbors with cost of sqrt(2); note adj8-adj4 are 8-neighbors that are not 4-neighbors

SID = sub2ind([h w], redgeCentX, redgeCentY); % starting point is the border point
FID = sub2ind([h w], row(indMinDist), col(indMinDist)); %end point is the closest point on bwMask

[cost,path] = DijkstraSparse(costMat, SID, FID);  
cutNodes = path; % last point (FID) is later dropped for value returned by function

cutNodes_up = cutNodes - 1; %nodes directly above cutNodes; not equivalent to the nodes_countercw defined later
cutNodes_down = cutNodes + 1; %nodes directly below cutNodes; not equivalent to the nodes_cw defined later

maskNodes = find(bwMask);
nodes_cw = setdiff(cutNodes_down, [FID; maskNodes]); %since we took shortcut to define nodes on cw and countercw, make sure FID and nodes on mask are not part of the paths
nodes_countercw = setdiff(cutNodes_up, [FID; maskNodes]);
cutNodes = setdiff(path, FID); %list of image nodes on the graph cut: excluding the closest point on the mask


% figure(88); subplot(1,2,1); imshow(bwMask); hold on; plot(redgeCentY, redgeCentX,'r*', col(indMinDist), row(indMinDist),' r*'); hold off
% subplot(1,2,2);  imshow(bwMask); hold on;
% [pathX pathY] = ind2sub([h w],cutNodes);
% [pathcwX pathcwY] = ind2sub([h w],nodes_cw);
% [pathcountercwX pathcountercwY] = ind2sub([h w],nodes_countercw);
% plot(pathY, pathX, 'r*', pathcwY, pathcwX,'y*', pathcountercwY,  pathcountercwX, 'b*');
% hold off;
% 
% figure(44);
% subplot(1,2,1); imshow(bwMask); hold on; plot(redgeCentY, redgeCentX,'r*', col(indMinDist), row(indMinDist),' r*'); hold off
% subplot(1,2,2); imshow(bwMask); hold on; 
% for i = 1: length(nodes_cw)
%         [cPointX cPointY] = ind2sub([h w], nodes_cw(i));
%         plot(cPointY, cPointX, 'y*');
%         
%         if i<=length(path)
%                 [cPointX cPointY] = ind2sub([h w], path(i));
%                 plot(cPointY, cPointX, 'r*');
%         end
% end
% hold off;

