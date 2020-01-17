function [individual] = initialize_kmeans( individual, numClusters, data )
%INITIALIZE_KMEANS Initialize one individual of the population using kmeans clustering
		[idx, clusters] = kmeans( data, numClusters, 'EmptyAction','singleton' ); 
		individual.bits(1:numClusters) = true;
		for k=1:numClusters
			individual.clusters(k).mean = clusters(k,:);
			idtmp = find( idx == k );
			individual.clusters(k).covariance = squareformSymmetric( cov([data(idtmp,:); data(idtmp,:)], 1) );
		end
end
