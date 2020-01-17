%Run GA-EM simulating the experiments described in Table 1 of Pernkopf05

clc
clear

addpath datasets
addpath GA-EM
addpath utils

load pendigits01234PCA.mat
load pendigitspca_lbl.mat

%omit unnecessary warnings
warning('off', 'stats:gmdistribution:FailedToConverge')
warning('off', 'stats:kmeans:FailedToConverge')


maxClusters=20;
pMutate=0.05;
pCross=0.8;
populationSize=6;
R=7;
tCorrel=0.95;
nRepetitions=10;
tAnnihil=10;

[numObjects, numFeatures] = size( pendigitsPCA );


disp('Starting EM with Random Initialization');
EM_Random.MDLs = inf .* ones(1, nRepetitions);
EM_Random.bestKs = zeros(1, nRepetitions);
EM_Random.EMSteps = zeros(1, nRepetitions);
EM_Random.Accs = zeros(1, nRepetitions);
EM_Random.confusionMatrix = [];
for r=1:nRepetitions
	disp(['Repetition ' num2str(r)])
	
	for k=2:maxClusters
		disp(['K: ' num2str(k)])
		objEM = gmdistribution.fit( pendigitsPCA, k, 'Options', statset('MaxIter', 10000, 'TolFun', 0.00001) );
		MDL = objEM.NlogL + ((k + k * numFeatures + k * numFeatures * (numFeatures+1)/2)/2 * log(numObjects));		
	
		if( MDL < EM_Random.MDLs(r) )
			EM_Random.MDLs(r) = MDL;
			EM_Random.bestKs(r) = objEM.NComponents;
			[ EM_Random.Accs(r) EM_Random.confusionMatrix(r).m ] = model_accuracy( posterior( objEM, pendigitsPCA ), pendigitsPCA_lbl);			 
		end
	
		EM_Random.EMSteps(r) = EM_Random.EMSteps(r) + objEM.Iters;
	end
end

disp('Starting GA-EM with Random Initialization');
GAEM_Random.MDLs = inf .* ones(1, nRepetitions);
GAEM_Random.bestKs = zeros(1, nRepetitions);
GAEM_Random.EMSteps = zeros(1, nRepetitions);
GAEM_Random.Accs = zeros(1, nRepetitions);
GAEM_Random.confusionMatrix = [];
for r=1:nRepetitions
	[ bestPartition EMStepsGAEM ] = GA_EM( pendigitsPCA, R , pMutate, pCross, populationSize, tCorrel, maxClusters, tAnnihil, 0);
		
	GAEM_Random.MDLs(r) = bestPartition.MDL;
	GAEM_Random.bestKs(r) = sum( bestPartition.bits );
	[ GAEM_Random.Accs(r) GAEM_Random.confusionMatrix(r).m ] = model_accuracy( evaluate_responsabilities( bestPartition, pendigitsPCA ), pendigitsPCA_lbl);

	GAEM_Random.EMSteps(r) = EMStepsGAEM;
	disp('.');
end


disp('Starting EM with Kmeans Initialization');
EM_Kmeans.MDLs = inf .* ones(1, nRepetitions);
EM_Kmeans.bestKs = zeros(1, nRepetitions);
EM_Kmeans.EMSteps = zeros(1, nRepetitions);
EM_Kmeans.Accs = zeros(1, nRepetitions);
EM_Kmeans.confusionMatrix = [];
for r=1:nRepetitions
	disp(['Repetition ' num2str(r)])
	
	for k=2:maxClusters
		disp(['K: ' num2str(k)])
		[pis means covs ] = EM_init_kmeans( pendigitsPCA', k );
		objEM = gmdistribution.fit( pendigitsPCA, k, 'Start', struct('mu', means', 'Sigma', covs, 'PComponents', pis),...
		   	'Options', statset('MaxIter', 10000, 'TolFun', 0.00001) );
		MDL = objEM.NlogL + ((k + k * numFeatures + k * numFeatures * (numFeatures+1)/2)/2 * log(numObjects));		
	
		if( MDL < EM_Kmeans.MDLs(r) )
			EM_Kmeans.MDLs(r) = MDL;
			EM_Kmeans.bestKs(r) = objEM.NComponents;
			[ EM_Kmeans.Accs(r) EM_Kmeans.confusionMatrix(r).m ] = model_accuracy( posterior( objEM, pendigitsPCA ), pendigitsPCA_lbl);			 
		end
	
		EM_Kmeans.EMSteps(r) = EM_Kmeans.EMSteps(r) + objEM.Iters;
	end
end


disp('Starting GA-EM with Kmeans Initialization');
GAEM_Kmeans.MDLs = inf .* ones(1, nRepetitions);
GAEM_Kmeans.bestKs = zeros(1, nRepetitions);
GAEM_Kmeans.EMSteps = zeros(1, nRepetitions);
GAEM_Kmeans.Accs = zeros(1, nRepetitions);
GAEM_Kmeans.confusionMatrix = [];
for r=1:nRepetitions
	[ bestPartition EMStepsGAEM ] = GA_EM( pendigitsPCA, R , pMutate, pCross, populationSize, tCorrel, maxClusters, tAnnihil, 1);
		
	GAEM_Kmeans.MDLs(r) = bestPartition.MDL;
	GAEM_Kmeans.bestKs(r) = sum( bestPartition.bits );
	[ GAEM_Kmeans.Accs(r) GAEM_Kmeans.confusionMatrix(r).m ] = model_accuracy( evaluate_responsabilities( bestPartition, pendigitsPCA ), pendigitsPCA_lbl);

	GAEM_Kmeans.EMSteps(r) = EMStepsGAEM;
	disp('.');
end

disp('Obtained results:')
disp('For EM with Random Initialization')
[minMDL idxMinMDL] = min( EM_Random.MDLs);
fprintf('Min MDL: %.0f\n', minMDL )
fprintf('Av MDL: %.0f (%.1f)\n', mean( EM_Random.MDLs ), std( EM_Random.MDLs ) )
fprintf('Best Comp: %d\n', EM_Random.bestKs( idxMinMDL ))
fprintf('Av Comp: %.1f\n', mean(EM_Random.bestKs) )
fprintf('EM-Steps: %.1f (%.1f)\n', mean(EM_Random.EMSteps), std(EM_Random.EMSteps) )
fprintf('%%: %.2f\n', 100*EM_Random.Accs( idxMinMDL ))

disp('---------------------------------------------')
disp('For GA-EM with Random Initialization')
[minMDL idxMinMDL] = min( GAEM_Random.MDLs);
fprintf('Min MDL: %.0f\n', minMDL )
fprintf('Av MDL: %.0f (%.1f)\n', mean( GAEM_Random.MDLs ), std( GAEM_Random.MDLs ) )
fprintf('Best Comp: %d\n', GAEM_Random.bestKs( idxMinMDL ))
fprintf('Av Comp: %.1f\n', mean(GAEM_Random.bestKs) )
fprintf('EM-Steps: %.1f (%.1f)\n', mean(GAEM_Random.EMSteps), std(GAEM_Random.EMSteps) )
fprintf('%%: %.2f\n', 100*GAEM_Random.Accs( idxMinMDL ))


disp('---------------------------------------------')
disp('For EM with k-means Initialization')
[minMDL idxMinMDL] = min( EM_Kmeans.MDLs);
fprintf('Min MDL: %.0f\n', minMDL )
fprintf('Av MDL: %.0f (%.1f)\n', mean( EM_Kmeans.MDLs ), std( EM_Kmeans.MDLs ) )
fprintf('Best Comp: %d\n', EM_Kmeans.bestKs( idxMinMDL ))
fprintf('Av Comp: %.1f\n', mean(EM_Kmeans.bestKs) )
fprintf('EM-Steps: %.1f (%.1f)\n', mean(EM_Kmeans.EMSteps), std(EM_Kmeans.EMSteps) )
fprintf('%%: %.2f\n', 100*EM_Kmeans.Accs( idxMinMDL ))

disp('---------------------------------------------')
disp('For GA-EM with k-means Initialization')
[minMDL idxMinMDL] = min( GAEM_Kmeans.MDLs);
fprintf('Min MDL: %.0f\n', minMDL )
fprintf('Av MDL: %.0f (%.1f)\n', mean( GAEM_Kmeans.MDLs ), std( GAEM_Kmeans.MDLs ) )
fprintf('Best Comp: %d\n', GAEM_Kmeans.bestKs( idxMinMDL ))
fprintf('Av Comp: %.1f\n', mean(GAEM_Kmeans.bestKs) )
fprintf('EM-Steps: %.1f (%.1f)\n', mean(GAEM_Kmeans.EMSteps), std(GAEM_Kmeans.EMSteps) )
fprintf('%%: %.2f\n', 100*GAEM_Kmeans.Accs( idxMinMDL ))

