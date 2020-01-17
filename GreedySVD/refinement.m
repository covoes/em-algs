function [individual] = refinement(individual, data)
%Refines a individual using EM

global EMSteps;
global DEBUG;

[numObjects numFeatures] = size(data);

statOpts = statset('TolFun', 1e-5);
%regularize value
regV = 1e-1;


%	info_individual(individual)
numClusters = individual.numClusters;

covs = zeros( numFeatures, numFeatures, numClusters );
means = individual.mean(1:numClusters,:);

%create gmdistribution 
for k=1:numClusters
	%adding regularization value to avoid ill conditioning
    covs(:,:,k) = squareformSymmetric( individual.covariance(k,:) ) + eye(numFeatures)*regV;
end
mixingCoefficients = individual.mixCoef(1:numClusters);


if DEBUG
	fprintf(DEBUG,'#REFINEMENT\nOLD INDIVIDUAL (%d):%s\n',i,info_individual(individual));
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

erro = 1;
while erro
	try
		objEM = gmdistribution.fit(data, numClusters, ...
			'Start', struct( 'mu', means, 'Sigma', covs, 'PComponents', mixingCoefficients ), ...
			'Options', statOpts, 'Regularize', regV);
		erro =0;
	catch err
		err.identifier
		err.message
		for k=1:numClusters
			%[a1,b1] = cholcov(covs(:,:,k))
			%adding regularization value to avoid ill conditioning
			%covs(:,:,k)
			covs(:,:,k) = covs(:,:,k) + eye(numFeatures)*regV;
		end
	end
end


EMSteps = EMSteps + objEM.Iters;

%update individual parameters
individual.mean(1:numClusters,:) = objEM.mu;
individual.mixCoef(1:numClusters) = objEM.PComponents;
individual.lastFitness = individual.fitness;
individual.fitness = fitnessFunc( 'mdl', objEM, numObjects );
for k=1:numClusters
    individual.covariance(k,:) = squareformSymmetric( objEM.Sigma(:,:,k) );
	individual.determinant(k) = det( objEM.Sigma(:,:,k) );
    %storing the squared mahalanobis distance
    dif = bsxfun(@minus, data, objEM.mu(k,:));
	individual.distance(:,k) = sum((dif / objEM.Sigma(:,:,k)) .* dif,2);
end


if DEBUG
	fprintf(DEBUG,'\nNEW INDIVIDUAL:%s\n', info_individual(individual));
end

