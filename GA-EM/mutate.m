function [ P ] = mutate( P, data, pMutation, maxClusters, a_min )
%MUTATE Perform uniform mutation over an individual
%
%	In the bit part of the individual just flip bits.
%	In the second part a random number is inserted in the means
%	of the individuals selected. These means are defined based on 
%	an upper and lower bound determined from the dataset (min/max values from data). 
%	The probability of modifying the second part is scale down by L.
%	The best solution is not modified.
%
%Parameters:
%P              : population of individuals
%data           : data NxM to be clustered
%pMutation      : probability of mutation
%maxClusters    : maximum number of clusters
%a_min	   	    : best solution found so far to avoid modifying it


	for i=1:length(P)

		%do not modify the best solution
		if( isequal(P(i).clusters, a_min.clusters) && isequal(P(i).bits, a_min.bits) )
			continue;
		end


		%define who is going to be mutated
		probs = rand(maxClusters, 1);
		idxMutation = find(probs <= pMutation);

		%flipping bit (part A)
		P(i).bits(idxMutation) = ~P(i).bits(idxMutation);

		%scale down pMutation by factor L
		numFeatures = size(data,2);
		L = numFeatures + numFeatures * (numFeatures + 1)/2;
		pMutation = pMutation/L;
		idxMutation = find(probs <= pMutation);

		%modify means (part B)
		mins = min(data);
		maxs = max(data);
		for k=1:length(idxMutation)
			P(i).clusters( idxMutation(k) ).mean = mins + (maxs - mins) .* rand();
		end

	end
	
end
