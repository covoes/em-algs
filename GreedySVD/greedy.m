function [bestPartition EMSteps tFinal] = greedy( data, maxClusters, maxIterations, maxKMSIter )
%GREEDY greedy algorithm based on PR03_Zhang
%
%Parameters:
%data              : dataset to be clustered NxM where N is the number of objects and M the number of features
%maxClusters       : maximum number of clusters to estimate
%maxIterations     : maximum number of iterations (termination criterion)
%maxKMSIter        : define the number of kmeans iteration used for initialization


global EMSteps;
global DEBUG;

%%Control Variables
%maximum candidates analyzed
max_candidates = 5; %value from ueda

%DEBUG must be a valid fid to be on and 0 to be off
DEBUG=0;
%DEBUG=fopen('logGREEDY.txt','w');
%EXTRA_INFO is used to save information during generations
%use 0 when information wont be needed and the name of .mat file 
%that will hold the info otherwise
EXTRA_INFO=0;
%EXTRA_INFO='extraInfo.mat';


EMSteps = 0;
lastBestFitness = 0;

if EXTRA_INFO
	INFO_K = NaN([maxIterations 1]);
	INFO_FIT = NaN([maxIterations 1]);
end

tIni = tic;

individual = initialize(data, maxClusters, maxKMSIter);

for t=1:maxIterations

	cur_fitness = individual.fitness;

	ranking_split = rank_split(individual, data);
	ranking_merge = rank_merge(individual, data);

%	fprintf('cur %g k%d\n',individual.fitness,individual.numClusters)
	modified=0;
	for c=1:max_candidates
		if individual.numClusters > 1 & c < length(ranking_merge)
			indiv_merge = merge(individual, data, ranking_merge, c);
			indiv_merge = refinement(indiv_merge, data);
%			fprintf('merg %g\n',indiv_merge.fitness)
			if (indiv_merge.fitness < cur_fitness)
				individual = indiv_merge;
				modified=1;
				break;
			end
		end

		%only test split if passed the merge (break)
		if individual.numClusters < maxClusters & c < individual.numClusters 			
			indiv_split = split(individual, data, ranking_split, c);
			indiv_split = refinement(indiv_split, data);
%			fprintf('split %g\n',indiv_split.fitness)
			if (indiv_split.fitness < cur_fitness)
				individual = indiv_split;
				modified=1;
				break;
			end
		end
	end
	%fprintf('-----------------------\n')

	if ~modified
		%no improvement was found
		break;
	end

	
	if EXTRA_INFO
		INFO_FIT(t,:) = [ individual.fitness ];
		INFO_K(t,:) = [ individual.numClusters ];
		save(EXTRA_INFO, 'INFO_FIT', 'INFO_K');
	end	
end

bestPartition = individual;

tFinal=toc(tIni);

if DEBUG
	fclose(DEBUG);
end


end
