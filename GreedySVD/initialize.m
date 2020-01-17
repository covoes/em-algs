function [individual] = initialize(data, maxClusters, maxKMSIter)
%Initialize the population with k-means

[numObjects numFeatures] = size(data);

individual = struct( 'mean', NaN([maxClusters numFeatures]), ...
					 'covariance', repmat(squareformSymmetric(NaN(numFeatures))', maxClusters, 1), ...
					 'mixCoef', NaN([1 maxClusters]), ...
					 'numClusters', NaN, ...
					 'fitness', NaN, ...
					 'lastFitness', NaN );


numClusters = randi([2 maxClusters], 1);
[idx, clusters] = kmeans( data, numClusters, 'EmptyAction', 'singleton','Options', statset('MaxIter',maxKMSIter)); 
individual.numClusters = numClusters;
for k=1:numClusters
	individual.mean(k,:) = clusters(k,:);
	idtmp = find( idx == k );
	if length(idtmp) > 0
		%the duplicate is because of the case nObjInCluster=1
		covariance = cov([data(idtmp,:); data(idtmp,:)]);
		individual.covariance(k,:) = squareformSymmetric( covariance );
		individual.mixCoef(k) = length(idtmp)/numObjects;
	else
		individual.mixCoef(k) = 0.1;
		individual.covariance(k) = squareformSymmetric(1e-5*eye(numFeatures));
	end
end
individual.mixCoef = individual.mixCoef./nansum(individual.mixCoef);

individual = refinement(individual, data);
