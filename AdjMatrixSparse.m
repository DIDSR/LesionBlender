function adj = AdjMatrixSparse(matSize, connectType)
% DIFFERENCE WITH ADJMATRIX IS THAT HERE THE OUTPUT OF THE FUNCTION IS A
% SPARSE MATRIX;
% function that returns a standard adjacency matrix for input matrix with
% size matSize, according to desired connection type (4 or 8); output is
% size prod(matSize)*prod(matSize)
% Taken from:
% http://stackoverflow.com/questions/3277541/construct-adjacency-matrix-in-matlab
r = matSize(2); %reason for switching described below
c = matSize(1);
% original logic used to  generate the adjacency matrix assumes row-wise node
%numbering:
% 1-2-3
% |   |   |
% 4-5-6
% |   |   |
% 7-8-9
% but we have column-wise numbered nodes, so height width indices are switched above 

% define a generic adjacency matrix based on 4 or 8 connectivity; size will
% be r*c by r*c
if connectType == 4 % logic for this part based on the regular patterns in adjacency matrix
        
        diagVec1 = sparse(repmat([ones(c-1,1); 0],r,1));  %# Make the first diagonal vector
        %#   (for horizontal connections)
        diagVec1 = diagVec1(1:end-1);             %# Remove the last value
        diagVec2 = sparse(ones(c*(r-1),1));               %# Make the second diagonal vector
        %#   (for vertical connections)
        adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
                diag(diagVec2,c);        
        adj = adj+adj.';                             %'# Add the matrix to a transposed
                                          %#   copy of itself to make it
                                          %#   symmetric
        
elseif connectType ==8
       
        diagVec1 = sparse(repmat([ones(c-1,1); 0],r,1));  %# Make the first diagonal vector
        %#   (for horizontal connections)
        diagVec1 = diagVec1(1:end-1);             %# Remove the last value
        diagVec2 = sparse([0; diagVec1(1:(c*(r-1)))]);    %# Make the second diagonal vector
        %#   (for anti-diagonal connections)
        diagVec3 = sparse(ones(c*(r-1),1));               %# Make the third diagonal vector
        %#   (for vertical connections)
        diagVec4 = sparse(diagVec2(2:end-1));             %# Make the fourth diagonal vector
        %#   (for diagonal connections)
        adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
                diag(diagVec2,c-1)+...
                diag(diagVec3,c)+...
                diag(diagVec4,c+1);
        adj = adj+adj.'; %'# Add the matrix to a transposed
                                          %#   copy of itself to make it
                                          %#   symmetric
end