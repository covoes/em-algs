function [ P ] = initialize( populationSize, data, maxClusters )
%INITIALIZE Initialize the structure of a set of individuals
%
%Parameters:
%populationSize : number of individuals in the population
%data           : data numObjects X numFeatures to be clustered
%maxClusters    : maximum number of clusters
[numObjects, numFeatures] = size( data );

cluster = struct( 'mean', zeros(1, numFeatures), ...
            'covariance', squareformSymmetric(eye(numFeatures)) );

P = struct( 'bits', 1:maxClusters > maxClusters, ...
            'clusters', repmat(cluster, maxClusters, 1), ...
            'fitness', inf);
P = repmat(P, populationSize, 1);

end
