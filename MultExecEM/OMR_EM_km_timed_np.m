function [bestGMMMin] = OMR_EM_km_timed_np(data, maxClusters, np)
%Ordered Multiple Runs of Expectation Maximization with a pre-specified number of partitions, and using MDL as model selection criteria

tIni = tic;
bestMDL = realmax;	
nObj = size(data,1);
EMSteps = 0;
bestGMMMin = {};
opts = statset('MaxIter',1000);
loop = 0;
optsk = statset('MaxIter',5);
regV=1e-5;
for n=1:np
	for k=2:maxClusters
		[idx, clusters] = kmeans( data, k, 'EmptyAction', 'singleton','Options', optsk );
		means = []; covariance = []; mixCoef = [];
		for kk=1:k
			means(kk,:) = clusters(kk,:);
			idtmp = find( idx == kk );
			if length(idtmp) > 0
				%the duplicate is because of the case nObjInCluster=1
				covariance(:,:,kk) = cov([data(idtmp,:); data(idtmp,:)]) + eye(size(data,2))*regV;
				mixCoef(kk) = length(idtmp)/size(data,1) + 0.1;
			else
				covariance(:,:,kk) = eye(size(data,2))*regV;
				mixCoef(kk) = 0.1;
			end
		end
		mixCoef = mixCoef/sum(mixCoef);
		initP = struct('mu', means, 'Sigma', covariance, 'PComponents', mixCoef);
		gmmObt = gmdistribution.fit(data,k,'Start', initP,'Regularize', regV,'Options',opts);
		mdlObt = fitnessFunc( 'mdl', gmmObt, nObj );
		EMSteps = EMSteps + gmmObt.Iters;
		if( mdlObt < bestMDL )
			bestGMM = gmmObt;
			bestMDL = mdlObt;
		end
	end
	tSeg = ceil(toc(tIni));
	bestGMMMin{n}.m = bestGMM;
	bestGMMMin{n}.tempoSeg = tSeg;
	loop = loop + 1;
end

end
