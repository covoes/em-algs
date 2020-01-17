function [P] = mutation(P, maxClusters, data)
%MUTATION Performs mutation on the indivuals in the population

global DEBUG;


for i=1:length(P)

	numClusters = P(i).numClusters;
	pElim = (numClusters-2)/(maxClusters-2);
	if (rand() > pElim)
		mutOpToApply = 'create';
	else
		mutOpToApply = 'elim';
	end

	if DEBUG
		fprintf(DEBUG,'\n\n\n-----%s\n\n#MUTATION\nOLD INDIVIDUAL (%d):%s\n',mutOpToApply,i,info_individual(P(i)));
	end

	
	[posterior gauss] = computePosterior( P(i), data );	
	oldIndiv = P(i);
	if strcmp(mutOpToApply,'elim')
		%remove clusters 
		lnTerms = bsxfun(@plus, log(eps+P(i).mixCoef(1:numClusters)), log(eps+gauss));
		%transform in positive so that higher values will have more probability to mutate
		valuesMutationClusters = -1*sum(posterior .* lnTerms, 1);
		[valuesSorted idxSorted] = sort( valuesMutationClusters );
		valuesMutationClusters( idxSorted ) = 1:numClusters;
		
		probs = valuesMutationClusters ./ sum(valuesMutationClusters);
		z = randi([1 numClusters-2]);
		%chosen are the clusters selected for removal
		chosen = roulette_without_reposition( probs, z );
		survivors = setdiff(1:numClusters, chosen);
		P(i).mean = P(i).mean(survivors,:);
		P(i).covariance = P(i).covariance(survivors,:);
		P(i).mixCoef = P(i).mixCoef(survivors);
		%re-normalize
		P(i).mixCoef = P(i).mixCoef ./ nansum(P(i).mixCoef);
		P(i).numClusters = length(survivors);
	else
		%create new clusters
		logsPost = log2(eps+posterior);
		objEntropies = -sum(posterior.*logsPost,2);
		
		[valuesSorted idxSorted] = sort( objEntropies );
		objEntropies( idxSorted ) = 1:length(objEntropies);

		probs = objEntropies./sum(objEntropies);
		z = randi([1 maxClusters-numClusters]);
		%chosen are the objects that will be used to create clusters
		chosen = roulette_without_reposition( probs, z );
		variances = var(data);
		P(i).numClusters = numClusters + z;
		%recover the clusters that were most probable to generate these points
		[~,oldClusters] = max(posterior(chosen,:), [], 2);
		for nc=1:z
			P(i).mean(numClusters+nc,:) = data(chosen(nc),:);
			P(i).covariance(numClusters+nc,:) = squareformSymmetric( 0.1*diag(variances) );
			P(i).mixCoef(numClusters+nc) = P(i).mixCoef(oldClusters(nc))/2;
			P(i).mixCoef(oldClusters(nc)) = P(i).mixCoef(oldClusters(nc))/2;
		end
	end

	if DEBUG
		assert(P(i).numClusters >= 2, 'Too few clusters')
		assert(P(i).numClusters <= maxClusters, 'Too many clusters')
		fprintf(DEBUG,'#MUTATION\nNEW INDIVIDUAL (%d):%s\n',i,info_individual(P(i)));
	end


end


