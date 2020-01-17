function [ responsabilities correspondentClusters ] = evaluate_responsabilities( individual, data ) 
%EVALUATE_RESPONSABILITIES Evaluates the responsabilities of each component
%
%Parameters:
%individual:  One individual of the population (a structure defined in "initialize")
%data      :  Dataset with numObjects X numFeatures dimensions

	numClusters = sum( individual.bits );
	correspondentClusters = find( individual.bits == 1);
	if(numClusters <= 1)
		%not much to do
		responsabilities = ones(size(data,1),1);
		return;
	end
	numFeatures = size( data, 2 );
	%evaluate posterior probabilities
	covs = zeros( numFeatures, numFeatures, numClusters );
	means = zeros( numClusters, numFeatures );
	for k=1:numClusters
		idxCluster = correspondentClusters( k );
		%1e-5 is used for regularization in case of singular covariance matrices
		covs(:,:,k) = squareformSymmetric( individual.clusters( idxCluster ).covariance ) + diag(repmat(1e-5, numFeatures,1));
		means(k,:) = individual.clusters( idxCluster ).mean;
	end
    mixingCoefficients = repmat(1/numClusters, numClusters, 1);

	if( isfield( individual, 'mixingCoefficients' ) )
		mixingCoefficients = individual.mixingCoefficients;
	end
    
	objEM = gmdistribution( means, covs, mixingCoefficients);
	responsabilities = posterior(objEM, data);
end
