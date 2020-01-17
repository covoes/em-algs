function [bestPartSegs] = FEAC_EM_timed_with_survivors( data, maxClusters, sizePopulation, ...
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


EMSteps = 0;

bestPartSegs = {};


P = initialize(data, maxClusters, sizePopulation,maxKMSIter);

while true

	P = refinement(P, data, maxEMIter, fitnessFName,maxKMSIter);

	%fprintf('nc: %.1f (%.1f)\n',mean([P(:).numClusters]),std([P(:).numClusters]));
	%fprintf('fit: %.1f (%.1f)\n',mean([P(:).fitness]),std([P(:).fitness]));


	[curBestFitness idx] = min([ P(:).fitness ]);
	bestPartition = P(idx);

	tSeg = 	ceil(toc(tIni));
	bestPartSegs{tSeg} = bestPartition;

	if(tSeg > maxTime)
		return;
	end


	Pmut = mutation(P, 0, splitFName, splitOpName, mergeFName, mergeOpName, maxKMSIter, maxClusters, data);
	Pmut = refinement(Pmut, data, maxEMIter, fitnessFName, maxKMSIter);
	P = selection_with_mutated_mulambda(P,Pmut,sizePopulation);

end


end
