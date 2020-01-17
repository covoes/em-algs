function [ P, a_min ] = perform_EM( P, data, R, a_min )
%perform_EM Perform R iterations of the EM algorithm in each individual of P and
%returns the coefficients of the best solution found 
%using data
%
%Parameters:
%P    : population
%data : numObjects X numFeatures  data matrix
%R    : number of EM iteration to execute
%a_min: Best individual so far
%
%Return:
%P    : updated population
%a_min: best individual found or the a_min given

global EMSteps

[numObjects, numFeatures] = size(data);

statOpts = statset('MaxIter', R, 'TolFun', 1e-5);

for i=1:length(P)

    numClusters = sum( P(i).bits );

	if( numClusters <= 1 )
		P(i).fitness = NaN;%realmax;
		continue;
	end

	covs = zeros( numFeatures, numFeatures, numClusters );
    means = zeros( numClusters, numFeatures );

	%an indexing is needed since it's possible that the cluster are not distributed contiguosly in 
	%the individual, correspondentClusters and idxCluster are used for that
	correspondentClusters = find( P(i).bits == 1 );
    for k=1:numClusters
		idxCluster = correspondentClusters(k);
        covs(:,:,k) = squareformSymmetric( P( i ).clusters( idxCluster ).covariance );
        means(k,:) = P( i ).clusters( idxCluster ).mean;
    end
    mixingCoefficients = repmat(1/numClusters, numClusters, 1);
    
	%if its refining the best partition the mixing coefficients are available
    if( isequal(a_min.clusters, P(i).clusters) && isequal(a_min.bits, P(i).bits ) )
        mixingCoefficients = a_min.mixingCoefficients;
    end

	
	objEM = gmdistribution.fit(data, numClusters, ...
        'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
        'Options', statOpts, 'Regularize', 1e-5);

	EMSteps = EMSteps + objEM.Iters;

    for k=1:numClusters
		idxCluster = correspondentClusters(k);
        P(i).clusters( idxCluster ).covariance = squareformSymmetric( objEM.Sigma(:,:,k) );
        P(i).clusters( idxCluster ).mean = objEM.mu(k,:);
    end

	%penalize only solutions with more than one cluster
	P(i).fitness = objEM.NlogL;
	if( numClusters > 1)
		%penalize using MDL criteria with relation to the free parameters
		penalMDL = (numClusters ... %mixing coefficients
			 + numClusters * numFeatures ... %means
		     + numClusters * numFeatures * (numFeatures+1)/2)/2 ... %covariance matrix
			 * log(numObjects); %data
		
		P(i).fitness = P(i).fitness + penalMDL;
	end

	if(P(i).fitness <= a_min.fitness)
		a_min = P(i);
		a_min.mixingCoefficients = objEM.PComponents;
		a_min.idv = i;
	end

end

end

