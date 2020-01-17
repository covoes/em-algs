function [P] = initialize(data, maxClusters, sizePopulation, maxKMSIter)
%Initialize the population with k-means

[numObjects numFeatures] = size(data);

individual = struct( 'mean', NaN([maxClusters numFeatures]), ...
					 'covariance', repmat(squareformSymmetric(NaN(numFeatures))', maxClusters, 1), ...
					 'mixCoef', NaN([1 maxClusters]), ...
					 'distance', NaN([numObjects maxClusters]), ...
					 'determinant', NaN([1 maxClusters]), ...
					 'numClusters', NaN, ...
					 'fitness', NaN);

P = repmat(individual, [1 sizePopulation]);

%initialize the parameter of each individual
opts = statset('MaxIter',maxKMSIter);
%changes from previous versions
%- now the individuals have the number of clusters linearly spaced between minimum and maximum number of clusters
if sizePopulation > 1
	numClusters = floor(linspace(2,maxClusters,sizePopulation));
else
	numClusters = randi([2 maxClusters],1);
end
for i=1:sizePopulation
	[idx, clusters] = kmeans( data, numClusters(i), 'EmptyAction', 'singleton','Options', opts ); 
	P(i).numClusters = numClusters(i);
	for k=1:numClusters(i)
		P(i).mean(k,:) = clusters(k,:);
		idtmp = find( idx == k );
		%the duplicate is because of the case nObjInCluster=1
		covariance = cov([data(idtmp,:); data(idtmp,:)]);
		P(i).covariance(k,:) = squareformSymmetric( covariance );
		P(i).mixCoef(k) = length(idtmp)/numObjects;
	end
end
