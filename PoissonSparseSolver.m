function dst = PoissonSparseSolver(src, Fh, Fv, msk)
tic;
lap = -grad2lap(Fh,Fv); %Fh & Fv are forward difference approximations, grad2lap retrieves the negative of laplacian based on forward and backward differences (see comments for grad2lap)
[row_msk col_msk] = find(msk); %if msk has more than two dimensions, it should have same value across the 3rd dimensions and taking "find" along one of them should be equivalent to any other 
minRowMsk = min(row_msk(:)); 
maxRowMsk = max(row_msk(:)); 
minColMsk = min(col_msk(:)); 
maxColMsk = max(col_msk(:)); 

if minRowMsk == 1 || maxRowMsk == size(src,1) || minColMsk == 1 || maxColMsk == size(src,2)
        error('PoissonSparseSolver: mask touches the sides of the image!');
end

src_crop = src(minRowMsk-1:maxRowMsk+1, minColMsk-1:maxColMsk+1);
lap_crop = lap(minRowMsk-1:maxRowMsk+1, minColMsk-1:maxColMsk+1);
msk_crop = msk(minRowMsk-1:maxRowMsk+1, minColMsk-1:maxColMsk+1);
[r c] = size(src_crop);
[r_knownValue c_knownValue] = find(~msk_crop); %find row/col pairs of pixels that aren't on the mask (i.e. pixels with known value)
indKnownValue = sub2ind([r c], r_knownValue, c_knownValue); %convert to linear index numbers
% r=800;c=800;

%initialize the A matrix in Ax=b for solving the poisson equation. Matrix A
%is same as 4-neighbor adjacency matrix, but with -4 along the diagonal for
%the pixels (nodes) located on msk_crop; After the matrix is initialized in
%this way, those rows corresponding to nodes not on the mask (i.e. nodes
%with known values) will be modified so that they get assigned their known
%values (they still need to be in A, because other nodes depend on their value).
r_fake = c; c_fake = r; %the algo for defining adjMatrix is based on row-wise numbering of nodes, whereas here we want col-wise numbering, so create new variables that switch row & columns (see comments in adjMatrix or AdjMatrixSparse)
diagVec0 = sparse(-2*ones(r_fake*c_fake,1)); %want the main diagonal to be equal to -4; since transpose later being added, set to -2 so that sum will be -4;

diagVec1 = sparse(repmat([ones(c_fake-1,1); 0],r_fake,1));  %# Make the first diagonal vector
%#   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);             %# Remove the last value
diagVec2 = sparse(ones(c_fake*(r_fake-1),1));               %# Make the second diagonal vector
%#   (for vertical connections)
adjs = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
        diag(diagVec2,c_fake) +...
        diag(diagVec0,0);
adjs = adjs+adjs.';

%initialize the b vector in Ax=b; the values for those rows in b that
%correspond to pixels with known values will be changed later. Will also
%change it into sparse representation later.
b = lap_crop(:);

%Now set the rows in A and b that correspond to pixels with known values
optionMode = 'col';
if strcmp(optionMode, 'row')
        for i = 1:length(indKnownValue)
                currentInd = indKnownValue(i);
                adjs(currentInd,:) = sparse(1, currentInd, 1, 1, r*c);
        end
elseif strcmp(optionMode, 'col') %accessing elements of sparse matrix in col format is much faster than row-wise access
        for i = 1:length(indKnownValue)
                currentInd = indKnownValue(i);
                adjs(:, currentInd) = sparse(currentInd, 1, 1, r*c, 1);
        end
        adjs = adjs'; %since col corresponding to an element was changed, take transpose to apply the change to the corresponding row for Ax=b;
end

b(indKnownValue) = src_crop(indKnownValue); %set values for pixels with known values
b = sparse(b);

dst_crop = adjs\b;
dst = src; %initialize
dst(minRowMsk-1:maxRowMsk+1, minColMsk-1:maxColMsk+1) = reshape(dst_crop, maxRowMsk+1 - minRowMsk+1 + 1, maxColMsk+1 -minColMsk+1 + 1);
disp(['Elapsed Time for Solution of Poisson Eqn via Sparse Solver: ' num2str(toc)]);


function lap = grad2lap(Fh, Fv)
% convert forward difference approximation to negative of laplacian
% approximation: the backward differences are Bh = circshift(Fh,[0 1]) &
% Bv=circshift(Fv,[1 0]), so laplacian kernel can be replicatd as sum of [1
% -2 1] and [1; -2 ; 1] convolved with image;
% the negative is returned because in the summation for Gauss-Siedel
% implementation here addition is used instead of subtraction
lap = (circshift(Fh,[0,1]) + circshift(Fv,[1,0]) - Fh - Fv);

