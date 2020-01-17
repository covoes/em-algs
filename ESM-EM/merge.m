function [individual]=merge(individual, mergeFName, mergeOpName, data)
%MERGE Merge two clusters in an individual
%
%mergeFName    : name of the function used to assess pairs of clusters to merge. Possible values are: 'cor', 'over'
%mergeOpName   : name of the operator that merge the clusters. Possible values are: 'simple'

global DEBUG;

%don't do nothing if there is only two clusters
if( individual.numClusters == 2)
	return
end

[numObjects numFeatures] = size(data);

switch( mergeFName )
	case 'cor'
		[posterior] = computePosterior( individual, data );
		%retrive the lower triangular part of the cosine similarity matrix in the vector form
		%is faster computing this way (even with the unnecessary multiplications) than with two fors
		%abs would not even be needed (since probabilities are nonnegative), 
		%but to be in agreement with the proposal...
		valuesMutationClusters = 1-abs(pdist(posterior','cosine'));

	case 'over'
		[posterior] = computePosterior( individual, data );
		[valsSorted idxClustSorted] = sort(posterior, 2, 'descend');
		overlap = zeros(individual.numClusters);
		%would be nice to make this code vectorized.. simple, apparently is not...
		for i=1:individual.numClusters
			for j=i+1:individual.numClusters
				%A = objects closer to Gaussian i with neighbor Gaussian j
				objsInA = find(idxClustSorted(:,i)==1 & idxClustSorted(:,j)==2);
				%B = objects closer to Gaussian j with neighbor Gaussian i
				objsInB = find(idxClustSorted(:,j)==1 & idxClustSorted(:,i)==2);

				overlap(i,j) = 0.5 * (sum( 1- (posterior(objsInA,i)- posterior(objsInA,j)))/length(objsInA) ...
					+ sum( 1- (posterior(objsInB,j)- posterior(objsInB,i)))/length(objsInB));
				overlap(j,i) = overlap(i,j);
			end
		end
		valuesMutationClusters = squareform(overlap);
end

%rank-based linear normalization
%assumes that higher values have more probability to mutate
[valuesSorted idxSorted] = sort( valuesMutationClusters );
pMutationClusters( idxSorted ) = 1:length(valuesMutationClusters);
pMutationClusters = pMutationClusters ./ sum(pMutationClusters);

%used for identification of the i,j-indexes of the pair of clusters
idxMat = squareform(1:length(idxSorted));

switch( mergeOpName )
	case 'simple'
		z = randi([1 floor(individual.numClusters/2)]);
		alreadyMutated = 1:individual.numClusters > inf;
		%used to identify the new indexes after the removal of some clusters
		oldNumClusters = individual.numClusters;
		convertIdx = 1:oldNumClusters;
		for c=1:z
			chosen = roulette( pMutationClusters, 1);
			%idxC1 and idxC2 are the indexes of the clusters that will be merged
			[idxC1 idxC2] = find( idxMat == chosen, 1 );
			
			%do not apply the operator twice in the same cluster 
			%it may be deleted or modified
			if (alreadyMutated(idxC1) || alreadyMutated(idxC2))
				continue;
			end
			alreadyMutated([ idxC1 idxC2 ]) = true;

			%convert old indexes to new ones
			idxC1 = convertIdx(idxC1);
			idxC2 = convertIdx(idxC2);

			%the new cluster will be saved in the smaller index
			%the greater index is removed and the other indexes are "decremented"
			idxNewCluster = idxC1;
			idxRemCluster = idxC2;
			if(idxC1 > idxC2)
				idxNewCluster = idxC2;
				idxRemCluster = idxC1;
			end

%			fprintf('merge de %d com %d - %d\n',idxNewCluster, idxRemCluster,c)
%			fprintf('\nantes')
%			info_individual(individual)

			meanC1 = individual.mean(idxC1,:);
			meanC2 = individual.mean(idxC2,:);
			mixCoefC1 = individual.mixCoef(idxC1);
			mixCoefC2 = individual.mixCoef(idxC2);
			%update the parameters of the "new" cluster
			individual.mixCoef(idxNewCluster) =  mixCoefC1 + mixCoefC2;
			individual.mean(idxNewCluster,:) = (mixCoefC1*meanC1 ...
						+ mixCoefC2*meanC2)/individual.mixCoef(idxNewCluster);
			covC1 = squareformSymmetric(individual.covariance(idxC1));
			covC2 = squareformSymmetric(individual.covariance(idxC2));
			tmpCov= (mixCoefC1 * covC1) + (mixCoefC1 * (meanC1' * meanC1))  ...
				  + (mixCoefC2 * covC2) + (mixCoefC2 * (meanC2' * meanC2));

			individual.covariance(idxNewCluster,:) = squareformSymmetric( (tmpCov ./ individual.mixCoef(idxNewCluster))...
			   	- (individual.mean(idxNewCluster,:)' * individual.mean(idxNewCluster,:)) );


			if DEBUG
				covs(:,:,1) = covC1;
				covs(:,:,2) = covC2;
				covL = squareformSymmetric(individual.covariance(idxNewCluster,:));
				meanL = individual.mean(idxNewCluster,:);
				%figure;
				%title('MERGE'); hold all;
				%plot(data(:,1),data(:,2),'.');
				%plotGMM([ meanC1' meanC2'], covs, [ .3 .3 .3 ], 1);
				%plotGMM( meanL' , covL, [ 0 .8 .0 ], 1);
			end


			%remove the other cluster and adjust the table for index conversion
			convertIdx(idxRemCluster+1:oldNumClusters) = convertIdx(idxRemCluster+1:oldNumClusters)-1;
			maxClusters = size(individual.mean,1);
			individual.mean(idxRemCluster:maxClusters,:) = [individual.mean(idxRemCluster+1:maxClusters,:); NaN([1 numFeatures])];
			individual.covariance(idxRemCluster:maxClusters,:) = [individual.covariance(idxRemCluster+1:maxClusters,:); NaN([1 length(individual.covariance(1,:))])];
			individual.mixCoef(idxRemCluster:maxClusters) = [individual.mixCoef(idxRemCluster+1:maxClusters)'; NaN];

			individual.numClusters=individual.numClusters-1;

%			fprintf('\ndps')
%			info_individual(individual)
		end
end
