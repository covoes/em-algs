function [bestGMM EMSteps] = OMR_EM(data, numberOfPartitions, maxClusters)
%Ordered Multiple Runs of Expectation Maximization with the outer loop for the partitions and the inner loop for the clusters. 
%Using MDL as model selection criteria

	bestMDL = realmax;	
	nObj = size(data,1);
	EMSteps = 0;
	for np=1:numberOfPartitions
		for k=2:maxClusters
			%tries to fit a gmm with randomly selected objects and variances along the dimensions of the initial covariance matrices
			gmmObt = gmdistribution.fit(data,k,'Regularize', 1e-5);
			mdlObt = fitnessFunc( 'mdl', gmmObt, nObj );
			EMSteps = EMSteps + gmmObt.Iters;
			if( mdlObt < bestMDL )
				bestGMM = gmmObt;
				bestMDL = mdlObt;
			end
		end
	end
end
