function [ P ] = initialize_initial_pop( populationSize, data, maxClusters, kmeansInitialization )
%INITIALIZE_INITIAL_POP Initialize the initial population using the procedure described in:
%Figueiredo & Jain 2002 or with values obtained with k-means
%   Initialize the population using the procedure described in: Figueiredo & Jain 
%Unsupervised Learning of Finite Mixture Models 2002. This is the strategy
%described as 1st strategy in Pernkop05. 
% The covariance are given as fraction of the global variance and the means are selected as random points of the data.
% The initialization by k-means is also implemented.
% The full individual is initialized to avoid problems in mutation.
%
%
%Parameters:
%populationSize      : number of individuals in the population
%data                : data numObjects X numFeatures to be clustered
%maxClusters         : maximum number of clusters
%kmeansInitialization: true if the individuals should be initialized with kmeans clustering

%each individual must have a different number of clusters so the
%initial population size may be bigger
P = initialize( max(populationSize, maxClusters-1), data, maxClusters );

[numObjects, numFeatures] = size( data );

numClusters = randsample([2:maxClusters], maxClusters-1);

%if the size of the population is bigger than the number of clusters
%than complete the rest of individual with a random number of clusters
if( populationSize > maxClusters-1 )
	numClusters( maxClusters:populationSize ) = 1 + randi( maxClusters-1, [1 populationSize-maxClusters+1] );
end


if kmeansInitialization

	for i=1:max(populationSize, maxClusters-1)	 		
		P(i) = initialize_kmeans( P(i), numClusters(i), data );
		%the rest of the means are initialized randomly to avoid problems in mutation
		for k=numClusters(i)+1:maxClusters
			P(i).clusters(k).mean = data( randi(numObjects) , :);			
		end
	end

else

	%initialization according to figueiredo
	globalMean = mean( data );
	%bsxfun is faster
	%diff = data - repmat(globalMean, numObjects, 1);
	diff = bsxfun(@minus, data, globalMean); 
	sigma_2 = (1/(10*numFeatures)) * trace( (diff' * diff)/numObjects );

	covMatrix = sigma_2 * eye(numFeatures);

	for i=1:max(populationSize, maxClusters-1)
		P(i).bits(1:numClusters(i)) = true;
		%initialize full individual
		seeds = randsample( numObjects, maxClusters);
		for k=1:maxClusters
			P(i).clusters(k).mean = data( seeds(k), :);
			P(i).clusters(k).covariance = squareformSymmetric(covMatrix);
		end
	end

end

end
