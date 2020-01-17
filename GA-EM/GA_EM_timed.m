function [ bestPartSegs ] = GA_EM_timed( data, R , pMutation, pRecombination,...
    populationSize, correlationThreshold, maxClusters, annihilationThreshold, kmeansInitialization, maxTime)
%GA_EM Implementation of the algorithm described in: Pernkopf & Bouchaffra: 2005
%   Algortihm described in Pernkopf & Bouchaffra: Genetic-Based EM Algorithm for Learning 
% Gaussian Mixture Models 2005
%
%%	 This version runs until maxTime (in minutes) is reached.
%%	 The results are returned in a struct similar to a population where the ith entry corresponds to 
%%	 the best individual found in i minutes
%
%
%Parameters:
%data                 : numObjectsXnumFeatures dataset
%R                    : number of EM steps to be performed in each solution
%pMutation            : probability of mutation operator 
%pRecombination       : probability of recombination operator
%populationSize       : number of individuals in the population
%correlationThreshold : maximum correlation accepted to not enforce mutation
%maxClusters          : maximum number of components in each individual 
%annihilationThreshold: minimun responsabilities sum that a component must have
%                       to not be annihilated
%kmeansInitialization : true if kmeans will be used for initializing individual,
%                       if false, figueiredo initialization will be used
%maxTime			  : maximum execution time in seconds

%number of EM steps performed
global EMSteps

tIni=tic();

EMSteps = 0;

%iterations counter
t = 0;
%number of components of the best solution 
OldSize = 0;
%number of generations without modification in the number of components
c_end = 0;
%initializing best solution values
a_min = initialize(1, data, maxClusters);

%initial population must have all the possible number of clusters
P = initialize_initial_pop( populationSize, data, maxClusters, kmeansInitialization);

bestPartSegs = {};

while true

    [P, a_min] = perform_EM( P, data, R, a_min );

	ncIdv=zeros(1,length(P));
	for idv=1:length(P)
		ncIdv(idv) = sum(P(idv).bits);
	end
	
	%DEBUG INF
	%fprintf('nc: %.1f (%.1f)\n',mean(ncIdv),std(ncIdv));
	%[minn minidx]=nanmin([P(:).fitness]);
	%fprintf('fit: %.3f (%.3f) :: %.3f (%d) -- %.3f (%d)\n',nanmean([P(:).fitness]./1200),nanstd([P(:).fitness]),P(a_min.idv).fitness/1200,a_min.idv, minn/1200,minidx );
	
	
	tSeg = ceil(toc(tIni));
	bestPartSegs{tSeg} = a_min;

	if(tSeg > maxTime)
		return;
	end

    childrens = recombine( P, data, maxClusters, pRecombination );
    
    [childrens, a_min] = perform_EM( childrens, data, R, a_min );


	P = selection( P, childrens, populationSize );


    P = enforce_mutation(P, data, correlationThreshold, a_min);

	P = mutate(P, data, pMutation, maxClusters, a_min);

    P = component_annihilation(P, data, annihilationThreshold, a_min);


	%fprintf('MDL min/avg/max %.2f/%.2f/%.2f a_min %.2f\n', min( [ P.MDL ] ),mean( [ P.MDL ] ),max( [ P.MDL ] ), a_min.MDL);
	%fprintf('K min/avg/max %.2f/%.2f/%.2f a_min %.2f\n', min( sum(vertcat(P.bits),2) ),mean( sum(vertcat(P.bits),2) ),max(sum(vertcat(P.bits),2) ), sum(a_min.bits));
	t = t + 1;
	
end

%run EM until the relative likelihood get smaller than 0.00001 (or 10 000 iterations)
[bestPartition unused] = perform_EM( a_min, data, 10000, a_min );

EMStepsRet = EMSteps;
end


