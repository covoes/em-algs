function [fitness] = fitnessFunc( FName, GMDistObj, numObjects, numClusters, numFeatures )
%Compute the fitness of a gmdistribution object.
%
%Possible values for FName are: 'mdl', 'fig','fig2'
%If adding a new function make sure that its in a minimization way
%i.e. the intention is to minimize the function

switch( lower(FName) )

	case 'mdl'
		penalMDL = (numClusters ... %mixing coefficients
			 + numClusters * numFeatures ... %means
		     + numClusters * numFeatures * (numFeatures+1)/2)/2 ... %covariance matrix
			 * log(numObjects); %data
		fitness = GMDistObj.NlogL + penalMDL;

	%Criterion defined in TPAMI02_Figueiredo 
	case 'fig'
		nparsover2 = (numFeatures + (numFeatures*(numFeatures+1)/2))/2;
		penalFIG = (nparsover2*(sum(log(GMDistObj.PComponents)))) + ((nparsover2 + 0.5)*numClusters*log(numObjects)) + (numClusters*nparsover2)+(numClusters/2);
		fitness = GMDistObj.NlogL + penalFIG;

	%Criterion defined in TPAMI02_Figueiredo according to its implementation
	case 'fig2'
		nparsover2 = (numFeatures + (numFeatures*(numFeatures+1)/2))/2;
		penalFIG = (nparsover2*(sum(log(GMDistObj.PComponents)))) + ((nparsover2 + 0.5)*numClusters*log(numObjects));
		fitness = GMDistObj.NlogL + penalFIG;


	otherwise
		error( 'Function %s not  implemented', FName );
	end
end
