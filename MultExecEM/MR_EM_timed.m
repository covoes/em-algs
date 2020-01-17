function [bestGMMMin] = MR_EM_timed(data, maxClusters, maxTime)
%Multiple Runs of Expectation Maximization using MDL as model selection criteria running for a specified number of minutes

tIni = tic;
bestMDL = realmax;	
nObj = size(data,1);
EMSteps = 0;
bestGMMMin = {};
opts = statset('MaxIter',10000000);
while true
	k = randi(maxClusters-1,1)+1;
	%tries to fit a gmm with randomly selected objects and variances along the dimensions of the initial covariance matrices
	gmmObt = gmdistribution.fit(data,k,'Regularize', 1e-5,'Options', opts );
	mdlObt = fitnessFunc( 'mdl', gmmObt, nObj );
	EMSteps = EMSteps + gmmObt.Iters;
	if( mdlObt < bestMDL )
		bestGMM = gmmObt;
		bestMDL = mdlObt;
	end
	tMin = ceil(toc(tIni)/60);
	bestGMMMin{tMin} = bestGMM;
	if(tMin > maxTime)
		return;
	end
end
end
