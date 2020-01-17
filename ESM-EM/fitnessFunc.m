function [fitness] = fitnessFunc( FName, GMDistObj, numObjects, numClusters, numFeatures )
%Compute the fitness of a gmdistribution object.
%
%Possible values for FName are: 'mdl'
%If adding a new function make sure that its in a minimization way
%i.e. the intention is to minimize the function

switch( lower(FName) )

	case 'mdl'
		penalMDL = (numClusters ... %mixing coefficients
			 + numClusters * numFeatures ... %means
		     + numClusters * numFeatures * (numFeatures+1)/2)/2 ... %covariance matrix
			 * log(numObjects); %data
		fitness = GMDistObj.NlogL + penalMDL;

	otherwise
		error( 'Function %s not  implemented', FName );
	end
end
