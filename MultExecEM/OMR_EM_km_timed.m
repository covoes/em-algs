function [bestGMMMin] = OMR_EM_km_timed(data, maxClusters, maxTime)
%Ordered Multiple Runs of Expectation Maximization with a pre-specified maximum allowed time in minutes, and using MDL as model selection criteria

tIni = tic;
bestMDL = realmax;	
nObj = size(data,1);
EMSteps = 0;
bestGMMMin = {};
opts = statset('MaxIter',10000000);
loop = 0;
optsk = statset('MaxIter',5);

while true
	for k=2:maxClusters
		[idx, clusters] = kmeans( data, k, 'EmptyAction', 'singleton','Options', optsk );
		means = []; covariance = []; mixCoef = [];
    	for kk=1:k
	        means(kk,:) = clusters(kk,:);
	        idtmp = find( idx == kk );
			if length(idtmp) > 0
		        %the duplicate is because of the case nObjInCluster=1
			    covariance(:,:,kk) = cov([data(idtmp,:); data(idtmp,:)]) + eye(size(data,2))*1e-5;
				mixCoef(kk) = length(idtmp)/size(data,1) + 0.1;
			else
			    covariance(:,:,kk) = eye(size(data,2))*1e-5;
				mixCoef(kk) = 0.1;
			end
		end
		mixCoef = mixCoef/sum(mixCoef);
		initP = struct('mu', means, 'Sigma', covariance, 'PComponents', mixCoef);
		%tries to fit a gmm with randomly selected objects and variances along the dimensions of the initial covariance matrices
		gmmObt = gmdistribution.fit(data,k,'Start', initP,'Regularize', 1e-5,'Options',opts);
		mdlObt = fitnessFunc( 'mdl', gmmObt, nObj );
		EMSteps = EMSteps + gmmObt.Iters;
		if( mdlObt < bestMDL )
			bestGMM = gmmObt;
			bestMDL = mdlObt;
		end
		tMin = ceil(toc(tIni)/60);
		bestGMMMin{tMin}.m = bestGMM;
		bestGMMMin{tMin}.loop = loop;
		if(tMin > maxTime)
			return;
		end
	end
	loop = loop + 1;
end

end
