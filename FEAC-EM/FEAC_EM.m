function [bestPartition EMSteps tFinal] = FEAC_EM( data, maxClusters, sizePopulation, ...
	   	maxGenerations, maxGenWOImprov, maxEMIter, fitnessFName, splitFName, ...
		mergeFName, splitOpName, mergeOpName, maxKMSIter )
%F-EAC algorithm adapted to use EM for partition refinement (not that simple)
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%maxClusters       : maximum number of clusters to estimate
%sizePopulation    : number of individuals in the population
%maxGenerations    : maximum number of generations (termination criterion)
%maxGenWOImprov    : maximum number of generations without improvement (termination criterion)
%maxEMIter         : maximum number of EM iterations for refinements in one EM execution
%fitnessFName      : name of fitness fuction to use. Possible are: 'mdl' (see fitnessFunc.m for details)
%splitFName        : name of the function used to assess clusters for splitting. Possible are: 'lkl','kul' (see split.m for details)
%mergeFName        : name of the function used to assess clusters for merging. Possible are: 'cor', 'over' (see merge.m for details)
%splitOpName       : name of the split operator to use. Possible are: 'svd', 'var' (see split.m for details)
%mergeOpName       : name of the merge operator to use. Possible are: 'simple' (see merge.m for details) 
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization


global EMSteps;
global DEBUG;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-4;
%number of individuals mantained by elitism
numElitists = 1;

%DEBUG must be a valid fid to be on and 0 to be off
DEBUG=0;
%DEBUG=fopen('logFEACEM.txt','w');
%EXTRA_INFO is used to save information during generations
%use 0 when information wont be needed and the name of .mat file 
%that will hold the info otherwise
EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';


EMSteps = 0;
genWOImprov = 0;
lastBestFitness = 0;


if EXTRA_INFO
	INFO_K = NaN([maxGenerations sizePopulation]);
	INFO_FIT = NaN([maxGenerations sizePopulation]);
end

tIni = tic;

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
			return
		end
	else
		genWOImprov = 0;
		lastBestFitness = curBestFitness;
	end

	P = selection(P, numElitists);
	P = mutation(P, numElitists, splitFName, splitOpName, mergeFName, mergeOpName, maxKMSIter, maxClusters, data);

%	g
%	bestPartition.fitness
	if EXTRA_INFO
		INFO_FIT(g,:) = [ P(:).fitness ];
		INFO_K(g,:) = [ P(:).numClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K');
	end	
end


tFinal=toc(tIni);

if DEBUG
	fclose(DEBUG);
end


end
