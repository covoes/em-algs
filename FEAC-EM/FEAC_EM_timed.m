function [bestPartMins] = FEAC_EM_timed( data, maxClusters, sizePopulation, ...
	   	maxEMIter, fitnessFName, splitFName, mergeFName, splitOpName, mergeOpName, maxKMSIter, maxTime )
%F-EAC algorithm adapted to use EM for partition refinement 
%
%
%%	 This version runs until maxTime (in minutes) is reached.
%%	 The results are returned in a struct similar to a population where the ith entry corresponds to 
%%	 the best individual found in i minutes
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%maxClusters       : maximum number of clusters to estimate
%sizePopulation    : number of individuals in the population
%maxEMIter         : maximum number of EM iterations for refinements in one EM execution
%fitnessFName      : name of fitness fuction to use. Possible are: 'mdl' (see fitnessFunc.m for details)
%splitFName        : name of the function used to assess clusters for splitting. Possible are: 'lkl','kul' (see split.m for details)
%mergeFName        : name of the function used to assess clusters for merging. Possible are: 'cor', 'over' (see merge.m for details)
%splitOpName       : name of the split operator to use. Possible are: 'svd', 'var' (see split.m for details)
%mergeOpName       : name of the merge operator to use. Possible are: 'simple' (see merge.m for details) 
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization
%maxTime


global EMSteps;
global DEBUG;


tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-4;
%number of individuals mantained by elitism
numElitists = 1;


EMSteps = 0;
genWOImprov = 0;
lastBestFitness = 0;

bestPartMins = {};


P = initialize(data, maxClusters, sizePopulation,maxKMSIter);

while true

	P = refinement(P, data, maxEMIter, fitnessFName,maxKMSIter);

	[curBestFitness idx] = min([ P(:).fitness ]);
	bestPartition = P(idx);

	tMin = 	ceil(toc(tIni)/60);
	bestPartMins{tMin} = bestPartition;

	if(tMin > maxTime)
		return;
	end

	%NOT USED IN TIMED VERSION
	%test termination criterion
%	if abs(lastBestFitness - curBestFitness) < tolerance
%		genWOImprov = genWOImprov + 1;
%		if genWOImprov == maxGenWOImprov
%			%no improvements in maxGenWOImprov iterations exit returning
			%the current best partition
%			tFinal=toc(tIni);
%			return
%		end
%	else
%		genWOImprov = 0;
%		lastBestFitness = curBestFitness;
%	end

	P = selection(P, numElitists);
	P = mutation(P, numElitists, splitFName, splitOpName, mergeFName, mergeOpName, maxKMSIter, maxClusters, data);

%	g
%	bestPartition.fitness
end


end
