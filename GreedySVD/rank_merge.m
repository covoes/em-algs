function [ranking]=rank_merge(individual, data)
%MERGE Merge two clusters in an individual
%
%mergeFName    : name of the function used to assess pairs of clusters to merge. Possible values are: 'cor', 'over'
%mergeOpName   : name of the operator that merge the clusters. Possible values are: 'simple'

global DEBUG;


[numObjects numFeatures] = size(data);

mergeFName='cor';

switch( mergeFName )
	case 'cor'
		[posterior] = computePosterior( individual, data );
		%retrive the lower triangular part of the cosine similarity matrix in the vector form
		%is faster computing this way (even with the unnecessary multiplications) than with two fors
		%abs would not even be needed (since probabilities are nonnegative), 
		%but to be in agreement with the proposal...
		%numClusters=size(posterior,2);
		%temp = zeros(numClusters);
		%for i=1:numClusters-1
		%	for j=i+1:numClusters
		%		tmp(i,j) = dot(posterior(:,i),posterior(:,j))/(norm(posterior(:,i))*norm(posterior(:,j)));
		%		tmp(j,i) = tmp(i,j);
		%	end
		%end
		%valuesMutationClusters = 1-abs(temp);
		valuesMutationClusters = 1-abs(pdist(posterior','cosine'));
end

%rank-based linear normalization
%assumes that higher values have more probability to mutate
[valuesSorted idxSorted] = sort( valuesMutationClusters );
ranking( idxSorted ) = length(valuesMutationClusters):-1:1;
end
