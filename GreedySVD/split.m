function [individual]=split(individual, data, ranking, candidate)
%SPLIT Splits a cluster in the individual.
%
%

global DEBUG;

%keeping the case switch, but not implementing other options, for now
splitOpName='svd';


%perform mutation using splitOpName
switch(splitOpName)
	case 'svd' %Singular Value Decomposition based Splitting

		chosen = ranking(candidate);
		%l is the cluster being splitted
		%it gives origin to clusters i and j 
		covL     = squareformSymmetric(individual.covariance( chosen, : ));
		meanL    = individual.mean( chosen,: );
		mixCoefL = individual.mixCoef( chosen );

		[U S V] = svd( covL );
		A = U*sqrt(S);
		a = A(:,1)';
		%these values are used in the experiments made in Zhang,2003
		alpha = 0.5;
		u = 0.5;
		betta = 0.5;

		%setting the new parameters
		newGroupI  = chosen;
		newGroupJ  = individual.numClusters+1;
		mixCoefI   = mixCoefL * alpha;
		mixCoefJ   = mixCoefL * (1-alpha);
		meanI      = meanL - sqrt( mixCoefJ/mixCoefI ) * u * a;
		meanJ      = meanL + sqrt( mixCoefI/mixCoefJ ) * u * a;
		aat        = a'*a;
		covI       = (mixCoefJ/mixCoefI) * covL + ( betta - betta * u^2 - 1 ) ...
					* (mixCoefL/mixCoefI) * aat + aat;			
		covJ       = (mixCoefI/mixCoefJ) * covL + ( betta * u^2 - betta - u^2 ) ...
					* (mixCoefL/mixCoefJ) * aat + aat;

		if DEBUG
			covs(:,:,1) = covI;
			covs(:,:,2) = covJ;
			%figure;
			%title('SPLIT'); hold all;
			%plot(data(:,1),data(:,2),'.'); 
			%plotGMM( meanL' , covL, [ .3 .3 .3 ], 1);
			%plotGMM([ meanI' meanJ'], covs, [ 0 .8 0 ], 1);
		end


		%updating individual
		individual.numClusters = individual.numClusters+1;
		individual.mean(newGroupI,:) = meanI;
		individual.mean(newGroupJ,:) = meanJ;
		individual.covariance(newGroupI,:) = squareformSymmetric( covI );
		individual.covariance(newGroupJ,:) = squareformSymmetric( covJ );
		individual.mixCoef(newGroupI) = mixCoefI;
		individual.mixCoef(newGroupJ) = mixCoefJ;

end

end
