function [ P ] = enforce_mutation( P, data, correlationThreshold,a_min )
%ENFORCE_MUTATION Modify some components if the correlation of their responsabilities with
% the responsabilities of other component is above the correlationThreshold
%   The modification is done either by removing the component or resetting it to
%   another data point, which strategy to use is selected randomly. The decision about 
%   which component (of the correlated pair) will be modified is random. 
%
%Parameters:
%P                   : population
%data                : numObjects X numFeatures dataset matrix
%correlationThreshold: maximun correlation allowed between the responsabilities of two gaussians
%a_min				 : best solution found so far to avoid modifying it

    numObjects = size(data,1);

	for i=1:length(P)

		%do not modify the best solution
		if( isequal(P(i).clusters, a_min.clusters) && isequal(P(i).bits, a_min.bits) )
			continue;
		end

		[ responsabilities correspondentClusters ] = evaluate_responsabilities( P(i), data );
        
        %find correlated pairs
        correlations = abs( corr(responsabilities) );
        tooCorrelated = triu( correlations > correlationThreshold , 1);
        
        [Is, Js] = find( tooCorrelated == 1 );

		candidateSet = zeros([length(Is) 1]);
		for p=1:length(Is)
			coin = rand;

            %choose cluster randomly
            candidateSet(p) = Is(p);
            if(coin > 0.5)
                candidateSet(p) = Js(p);
            end

		end

		%remove duplicates
		candidateSet = unique( candidateSet );

        %for each cluster in the candidateSet make a decision
        for p=1:length(candidateSet)
            coin = randi(2);
			idxCluster = correspondentClusters( candidateSet(p) );

			%make modification
            if(coin == 1)
                %remove component
                P(i).bits( idxCluster ) = 0;
            else
                %assign mean to a random object
				P(i).clusters( idxCluster ).mean = data( randi(numObjects), :);
            end
        end

	end
end

