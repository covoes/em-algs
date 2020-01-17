function [ranking]=rank_split(individual, data)
%RANK_SPLIT Generates a ranking of cluster to be splitted
%
%

%compute mutation probabilities for each cluster using splitFName
%assumes that higher values have more probability to mutate

global DEBUG;

%keeping the case switch, but not implementing other options, for now
splitFName='kul';

switch(splitFName)
	case 'kul'
		[posterior gauss] = computePosterior( individual, data );
		%assuming no objects in the same position.. for continuous features this is expected
		localDensity = bsxfun(@rdivide,posterior,sum(posterior, 1));
		%make sure it sums up to 1
		gauss = bsxfun(@rdivide, gauss, sum(gauss,1));
		valuesMutationClusters = localDensity .* log2( localDensity ./ gauss );
		valuesMutationClusters = sum(valuesMutationClusters, 1);		
end

%rank-based linear normalization
[valuesSorted idxSorted] = sort( valuesMutationClusters );
ranking( idxSorted ) = individual.numClusters:-1:1;
end
