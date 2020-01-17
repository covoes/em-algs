function [obj] = resultsGAEMToGMM(resultsStruct)
	nc = sum(resultsStruct.bits);
	nf = size(resultsStruct.clusters(1).mean,2);
	Me = zeros(nc,nf);
	Co = zeros(nf,nf,nc);
	idxC=0;
	%check if the NaN were not removed (in this case, the bits part was not updated)
	if(~any(isnan(resultsStruct.clusters(1).mean)))
		resultsStruct.bits = ones(1, length(resultsStruct.clusters));
	end
	for i=1:length(resultsStruct.bits)
		if(resultsStruct.bits(i) == 1)
			idxC = idxC + 1;
			Me(idxC,:) = resultsStruct.clusters(i).mean;
			Co(:,:,idxC) = squareformSymmetric(resultsStruct.clusters(i).covariance);
		end
	end
	W = resultsStruct.mixingCoefficients(1:nc);
	obj = gmdistribution(Me,Co,W);
end
