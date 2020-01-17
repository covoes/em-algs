function [bestPartition EMSteps tFinal] = ESMEM( data, maxClusters, sizePopulation, ...
	   	maxGenerations, maxGenWOImprov, maxEMIter, fitnessFName, splitFName, ...
		mergeFName, splitOpName, mergeOpName, maxKMSIter )
% Evolutionary Split & Merge Expectation Maximization
% Reference:
% Evolving Gaussian Mixture Models with Splitting and Merging Mutation Operators
%   Thiago Ferreira Covões, Eduardo Raul Hruschka, and Joydeep Ghosh
%   Evolutionary Computation 2016 24:2, 293-317
%
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%maxClusters       : maximum number of clusters to estimate
%sizePopulation    : number of individuals in the population
%maxGenerations    : maximum number of generations
%maxGenWOImprov    : maximum number of generations without improvement
%maxEMIter         : maximum number of EM iterations for refinements in one EM execution
%fitnessFName      : name of fitness fuction to use. Possible are: 'mdl' (see fitnessFunc.m for details)
%splitFName        : name of the function used to assess clusters for splitting. Possible are: 'lkl','kul' (see split.m for details)
%mergeFName        : name of the function used to assess clusters for merging. Possible are: 'cor', 'over' (see merge.m for details)
%splitOpName       : name of the split operator to use. Possible are: 'svd', 'var' (see split.m for details)
%mergeOpName       : name of the merge operator to use. Possible are: 'simple' (see merge.m for details)
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization


global EMSteps;
global DEBUG;

EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';

tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-5;


EMSteps = 0;
genWOImprov = 0;
lastBestFitness = 0;



P = initialize(data, maxClusters, sizePopulation,maxKMSIter);
for g=1:maxGenerations
	P = refinement(P, data, maxEMIter, fitnessFName,maxKMSIter);

	[curBestFitness idx] = min([ P(:).fitness ]);

	bestPartition = P(idx);

	%test termination criterion
	if abs(lastBestFitness - curBestFitness) < tolerance
		genWOImprov = genWOImprov + 1;
		if genWOImprov == maxGenWOImprov
			%no improvements in maxGenWOImprov iterations exit returning
   		%the current best partition
			tFinal=toc(tIni);
			g
			return
		end
	else
		genWOImprov = 0;
		lastBestFitness = curBestFitness;
	end

	Pmut = mutation(P, 0, splitFName, splitOpName, mergeFName, mergeOpName, maxKMSIter, maxClusters, data);
	Pmut = refinement(Pmut, data, maxEMIter, fitnessFName, maxKMSIter);
	P = selection_with_mutated_mulambda(P,Pmut,sizePopulation);
    %fprintf('\n-%.2f (%.2f)\n',mean([ P(:).numClusters ]), std([P(:).numClusters]))

%	g
%	bestPartition.fitness
	if EXTRA_INFO
		INFO_FIT(g,:) = [ P(:).fitness ];
		INFO_K(g,:) = [ P(:).numClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K');
	end

end

tFinal = toc(tIni);

end
