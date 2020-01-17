function [posterior gauss] = computePosterior( individual, data )
%COMPUTEPOSTERIOR Computer the posterior probability for each individual to each cluster
%
%In the second argument the gaussian pdfs for the objects are returned in
%a numObjects x numClusters matrix

[numObjects numFeatures] = size(data);
numClusters = individual.numClusters;
mixCoef = individual.mixCoef(1:numClusters);

posterior = NaN([numObjects numClusters]);

gauss = NaN([numObjects numClusters]);

const = (2*pi)^(numFeatures/2);
for k=1:numClusters
	gauss(:,k) = (1/(const* (individual.determinant(k))^0.5)) .* exp(-0.5.*individual.distance(:,k));
	posterior(:,k) = mixCoef(k) .* gauss(:,k);
end

%normalizing posteriors
sums = sum(posterior, 2);
posterior = bsxfun(@rdivide, posterior, sums);

end
