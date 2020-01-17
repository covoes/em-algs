function [ accuracy confusionMatrix ] = model_accuracy( responsabilities, labels )
%MODEL_ACCURACY Compute the accuracy of a clustering model according 
%to a pre-defined hard-clustering
%
% Each label in labels can be represented by a different component (cluster) 
% represented in responsabilities
%Parameters:
%responsabilities: numObjects X numClusters matrix of posterior probabilities
%labels:           numObjects X 1 column vector with the correct label for each
%                  object
%
%Return:
%accuracy
%confusionMatrix: numClusters X numLabels where each cell (i,j) is the number of objects
%                clustered to the ith cluster and labeled to the jth label

[numObjects, numClusters] = size(responsabilities);

numLabels = length(unique(labels));

confusionMatrix = zeros( numClusters, numLabels );

[junk clusters] = max( responsabilities, [], 2 );

for i=1:numClusters
	for j=1:numLabels
		confusionMatrix(i,j) = sum( clusters == i & labels == j );
	end
end

%assume that each cluster correspond to the label with most objects in the cluster
%allowing that more than one cluster represent the same label
numCorrect = max( confusionMatrix, [], 2);
accuracy =  sum( numCorrect ) / numObjects;


end
