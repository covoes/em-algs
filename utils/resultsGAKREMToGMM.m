function [obj] = resultsGAKREMToGMM(resultsStruct, data)
	[Co Me W] = estimate_parameters_with_kmeans(resultsStruct.medoids, data);
	obj = gmdistribution(Me,Co,W);
end
