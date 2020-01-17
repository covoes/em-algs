function [bestPartition EMSteps tFinal g] = ECE_EM( data, maxClusters, sizePopulation, ...
	   	maxGenerations, maxGenWOImprov, maxEMIter, fitnessFName, maxKMSIter )
%ECE_EM algorithm, Evolutionary Create & Eliminate EM  
%
%
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%maxClusters       : maximum number of clusters to estimate
%sizePopulation    : number of individuals in the population
%maxGenerations    : maximum number of gerenerations
%maxGenWOImprov    : maximum number of generations without improvement
%maxEMIter         : maximum number of EM iterations for refinements in one EM execution
%fitnessFName      : name of fitness fuction to use. Possible are: 'mdl' (see fitnessFunc.m for details)
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization


global EMSteps;
global DEBUG;

EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';


DEBUG=0;
%DEBUG=fopen('debug.txt','w');

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
	Pmut = mutation(P, maxClusters, data);
	Pmut = refinement(Pmut, data, maxEMIter, fitnessFName, maxKMSIter);
	P = selection_with_mutated_mulambda(P,Pmut,sizePopulation);
	%fprintf('\n-%.2f (%.2f)\n',mean([ P(:).numClusters ]), std([P(:).numClusters]))
	if EXTRA_INFO
		INFO_FIT(g,:) = [ P(:).fitness ];
		INFO_K(g,:) = [ P(:).numClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K');
	end	

end

tFinal = toc(tIni);

end
