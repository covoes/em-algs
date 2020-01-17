function [individual]=merge(individual, data, ranking, candidate)
%MERGE Merge two clusters in an individual
%

global DEBUG;

%don't do nothing if there is only two clusters
if( individual.numClusters == 2)
	return
end

[numObjects numFeatures] = size(data);

%used for identification of the i,j-indexes of the pair of clusters
idxMat = squareform(1:length(ranking));

%keeping the switch case for future use
mergeOpName='simple';

switch( mergeOpName )
	case 'simple'
		%idxC1 and idxC2 are the indexes of the clusters that will be merged
		[idxC1 idxC2] = find( idxMat == ranking(candidate), 1 );
			

		%the new cluster will be saved in the smaller index
		%the greater index is removed and the other indexes are "decremented"
		idxNewCluster = idxC1;
		idxRemCluster = idxC2;
		if(idxC1 > idxC2)
			idxNewCluster = idxC2;
			idxRemCluster = idxC1;
		end


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


		%remove the other cluster
		maxClusters = size(individual.mean,1);
		individual.mean(idxRemCluster:maxClusters,:) = [individual.mean(idxRemCluster+1:maxClusters,:); NaN([1 numFeatures])];
		individual.covariance(idxRemCluster:maxClusters,:) = [individual.covariance(idxRemCluster+1:maxClusters,:); NaN([1 length(individual.covariance(1,:))])];
		individual.mixCoef(idxRemCluster:maxClusters) = [individual.mixCoef(idxRemCluster+1:maxClusters)'; NaN];

		individual.numClusters=individual.numClusters-1;

end


end
