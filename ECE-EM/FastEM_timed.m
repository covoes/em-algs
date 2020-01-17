function [bestPartMins] = FastEM_timed( data, maxClusters, sizePopulation, ...
	   	maxEMIter, fitnessFName, maxKMSIter, maxTime )
%FASTEM algorithm  
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
%maxKMSIter        : if the variance based splitting is used, define the number of kmeans iteration used for parameter estimation, also used in the initialization
%maxTime


global EMSteps;
global DEBUG;


tIni = tic;

%%Control Variables
%tolerance to consider no improvement in EM
tolerance = 1e-4;


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


	Pmut = mutation(P, maxClusters, data);
	Pmut = refinement(Pmut, data, maxEMIter, fitnessFName, maxKMSIter);
	P = selection_with_mutated_mulambda(P,Pmut,sizePopulation);

end


end
