function [ Y ] = squareformSymmetric( X )
%SQUAREFORMSYMMETRIC Function similar to the squareform function but works for
%any symmetric matrix

if( ~isvector(X) )%(size(X, 2) > 1) 
    %transform matrix into a vector
    %copy only positions in the lower triangle
    idx = find(tril(ones(size(X))));    
    Y = X(idx);   
else
	X = X(:);
    %transform vector into matrix
    %X has x*(x+1)/2 elements where x is the order of the original matrix, the
    %original order can be computed as:
    Y = zeros( (-1 + sqrt( 1 + 8 * size(X,1) ))/2);
    idx = find(tril(ones(size(Y))));
    Y(idx) = X;
    Y = Y + Y' - diag(diag(Y));
end


end

