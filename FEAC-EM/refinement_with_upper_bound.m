function [P] = refinement_with_upper_bound(P, data, maxEMIter, fitnessFName)
%Refines each individual of the population using EM, however, uses a theoretical upper bound (ICML08_Zhang) to avoid unnecessary EM iterations on bad solutions, for this, at each EM iteration we see if the upper bound is less than the best fitness found, if it is, we stop iterating, otherwise we iterate until the llk difference be smaller than 1e-4

global EMSteps;
global DEBUG;

[numObjects numFeatures] = size(data);

statOpts = statset('MaxIter', 1);
%regularize value
regV = 1e-2;

curBestFitness = min([ P(:).fitness ]);

for i=1:length(P)


%	info_individual(P(i))
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
		%f=figure;
		%title(sprintf('REF %s',mat2str(mixingCoefficients))); hold all;
		%plot(data(:,1),data(:,2),'.'); 
		%plotGMM( means' , covs, [ .3 .3 .3 ], 1);
		%hold off;
		%[tmp d] =system('date +%s');
		%print(f, '-depsc2', sprintf('figs/ref_%d_%d.eps', i, mod(str2num(d),100000)));
		%close(f);
	end
	oldLLK = NaN;
	while true
		objEM = gmdistribution.fit(data, numClusters, ...
		    'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
			'Options', statOpts, 'Regularize', regV);
		EMSteps = EMSteps + objEM.Iters;

		fitness = fitnessFunc( fitnessFName, objEM, numObjects );
		%add upper bound at this point (see ICML08_Zhang) guardar likelihood por k? 
		delta = min([1 sqrt( (6 *(curBestFitness-(-objEM.NlogL)))/(n*min(objEM.PComponents)))]);
		upperBoundFitness = fitness + numObject * min(objEM.PComponents) * (delta^2)/6;
		%if we cant get better exit, otherwise continue iterating
		if upperBoundFitness > curBestFitness
			break;
		else
			%did not change enough, stop anyway
			if ((~isnan(oldLLK)) & (abs(oldLLK-objEM.NlogL) < 1e-4))
				break;
			else
				means = objEM.mu;
				mixingCoefficients = objEM.PComponents;
				covs = objEM.Sigma;
				oldLLK = objEM.NlogL;
			end
		end
	end
	%update individual parameters
	P(i).mean(1:numClusters,:) = objEM.mu;
	P(i).mixCoef(1:numClusters) = objEM.PComponents;
	P(i).lastFitness = P(i).fitness;
	P(i).fitness = fitnessFunc( fitnessFName, objEM, numObjects );
	for k=1:numClusters
        P(i).covariance(k,:) = squareformSymmetric( objEM.Sigma(:,:,k) );
		P(i).determinant(k) = det( objEM.Sigma(:,:,k) );
		%storing the squared mahalanobis distance
		dif = bsxfun(@minus, data, objEM.mu(k,:));
		P(i).distance(:,k) = sum((dif / objEM.Sigma(:,:,k)) .* dif,2);		
	end


	if DEBUG
		fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(P(i)));
	end

end

end
