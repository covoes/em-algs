function [obj] = resultsToGMM(resultsStruct)
	nc = resultsStruct.numClusters;
	nf = size(resultsStruct.mean,2);
	Me = resultsStruct.mean(1:nc,:);
	Co = zeros(nf,nf,nc);
	for i=1:nc
		Co(:,:,i) = squareformSymmetric(resultsStruct.covariance(i,:));
	end
	W = resultsStruct.mixCoef(1:nc);
	obj = gmdistribution(Me,Co,W);
end
