function [ P ] = component_annihilation( P, data, annihilationThreshold, a_min )
%COMPONENT_ANNIHILATION Eliminates components that have the responsabilities sum smaller than annihilationThreshold
%
%Parameters:
%P					  : population
%data                 : numObjectsXnumFeatures dataset
%annihilationThreshold: minimun responsabilities sum that a component must have
%                       to not be annihilated
%a_min				  : best solution found so far to avoid modifying it


	for i=1:length(P)

		%do not modify the best solution
		if( isequal(P(i).clusters, a_min.clusters) && isequal(P(i).bits, a_min.bits) )
			continue;
		end

		[ responsabilities componentsCorrespondent ] = evaluate_responsabilities( P(i), data );
		sums = sum(responsabilities, 1);
		annihilateSet = sums < annihilationThreshold;
		P(i).bits( componentsCorrespondent( annihilateSet ) ) = 0;
	end
end

