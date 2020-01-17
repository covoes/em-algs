function [bestGMMMin] = OMR_EM_timed(data, maxClusters, maxTime)
%Ordered Multiple Runs of Expectation Maximization with a pre-specified maximum allowed time in minutes, and using MDL as model selection criteria

tIni = tic;
bestMDL = realmax;	
nObj = size(data,1);
EMSteps = 0;
bestGMMMin = {};
opts = statset('MaxIter',10000000);
while true
	for k=2:maxClusters
		%tries to fit a gmm with randomly selected objects and variances along the dimensions of the initial covariance matrices
		gmmObt = gmdistribution.fit(data,k,'Regularize', 1e-5,'Options',opts);
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

end
