function [P] = selection(P, numElitists)
%SELECTION Performs elitism and roulette wheel selection on the population
%
%Assumes that the intent is to minimize the fitness function

global DEBUG;

	fitness = [P(:).fitness];
	sizePopulation = length(P);

	%elitism
	[bestFit idxSorted] = sort(fitness,'ascend');
	best = idxSorted(1:numElitists);
	
	%roulette
	ranks = zeros(sizePopulation,1);
	ranks(idxSorted) = sizePopulation:-1:1;
	chance = ranks/sum(ranks);	
	chosen = roulette( chance, sizePopulation-numElitists );
	
	if DEBUG
		fprintf(DEBUG,'#SELECTION\nFITNESS:\n\t%s',mat2str(fitness,4));
		fprintf(DEBUG,'\nNUM K:\n\t%s',mat2str([P(:).numClusters],4));
		fprintf(DEBUG,'\nPROBS:\n\t%s',mat2str(chance,4));
		fprintf(DEBUG,'\nELITISM:\n\t%s',mat2str(best));
		fprintf(DEBUG,'\nCHOSEN:\n\t%s\n',mat2str(chosen,4));
	end

	P(1:numElitists) = P(best);
	P(numElitists+1:sizePopulation) = P( chosen );

end
