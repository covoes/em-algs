function [P] = refinement(P, data, maxEMIter, fitnessFName, maxKMSIter)
%Refines each individual of the population using EM

global EMSteps;
global DEBUG;

[numObjects numFeatures] = size(data);

statOpts = statset('MaxIter', maxEMIter, 'TolFun', 1e-5);
%regularize value
regV = 1e-5;

for i=1:length(P)


	numClusters = P(i).numClusters;
		
	covs = zeros( numFeatures, numFeatures, numClusters );
    means = P(i).mean(1:numClusters,:);

	%create gmdistribution 
    for k=1:numClusters
		%adding regularization value to avoid ill conditioning
        covs(:,:,k) = squareformSymmetric( P( i ).covariance(k,:) ) + eye(numFeatures)*regV;
    end
    mixingCoefficients = P(i).mixCoef(1:numClusters);


	if DEBUG
		fprintf(DEBUG,'#REFINEMENT\nOLD INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
		save temp.mat
	end


	%if we find a problem with ill conditioned covariance matrices, restart a new solution
	try
		objEM = gmdistribution.fit(data, numClusters, ...
		    'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
			'Options', statOpts, 'Regularize', regV);
	catch err		
		fprintf('--%s\n',err.identifier)
		P(i) = initialize(data, length(P(i).mixCoef), 1, maxKMSIter);
		numClusters = P(i).numClusters;
		covs = zeros( numFeatures, numFeatures, numClusters );
	    means = P(i).mean(1:numClusters,:);

		%create gmdistribution 
	    for k=1:numClusters
			%adding regularization value to avoid ill conditioning
		    covs(:,:,k) = squareformSymmetric( P( i ).covariance(k,:) ) + eye(numFeatures)*regV;
	    end
		mixingCoefficients = P(i).mixCoef(1:numClusters);
		
		objEM = gmdistribution.fit(data, numClusters, ...
		    'Start', struct( 'mu', means, 'Sigma', covs , 'PComponents', mixingCoefficients ), ...
			'Options', statOpts, 'Regularize', regV);
	end

	EMSteps = EMSteps + objEM.Iters;
	
	%update individual parameters
	P(i).mean(1:numClusters,:) = objEM.mu;
	P(i).mixCoef(1:numClusters) = objEM.PComponents;
	P(i).fitness = fitnessFunc( fitnessFName, objEM, numObjects, numClusters, numFeatures );
	covs = objEM.Sigma;
	for k=1:numClusters
        P(i).covariance(k,:) = squareformSymmetric( covs(:,:,k) );
		P(i).determinant(k) = det( covs(:,:,k) );
		%storing the squared mahalanobis distance
		dif = bsxfun(@minus, data, P(i).mean(k,:));
		P(i).distance(:,k) = sum((dif / covs(:,:,k)) .* dif,2);		
	end

	if DEBUG
		fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(P(i)));
	end

end

end
